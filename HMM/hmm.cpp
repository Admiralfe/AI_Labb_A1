#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <algorithm>

#include "hmm.h"
#include "globals.h"
#include "matrix.h"

#define NO_OBS 9
#define NO_HS 3
#define VARIANCE 0.04
#define TIME_OUT -2

//define this constant to remove all asserts
#define SAFETY_OFF_HMM

using namespace globals;
using namespace std;

typedef matrix mat;
typedef vector<number> vec;

//create a new HMM with nearly uniform A, B and pi
Lambda::Lambda() {
    A = matrix::random_uniform(NO_HS, NO_HS, VARIANCE);
    B = matrix::random_uniform(NO_HS, NO_OBS, VARIANCE);
    pi = matrix::random_uniform(1, NO_HS, VARIANCE).get_row(0);
}

Lambda::Lambda(const matrix& transition, const matrix& emission, const vector<number>& init_state) {
    A = transition;
    B = emission;
    pi = init_state;
}

//Resets the model parameters of lambda to random_uniform values but keeps the observation sequence.
void Lambda::reset() {
    A = matrix::random_uniform(NO_HS, NO_HS, VARIANCE);
    B = matrix::random_uniform(NO_HS, NO_OBS, VARIANCE);
    pi = matrix::random_uniform(1, NO_HS, VARIANCE).get_row(0);
}

//guess the next observation based on a known sequence and an HMM
//also puts the estimated probability in max_prob
int hmm::next_obs_guess(const Lambda& lambda, const pair<vector<int>, int>& observations, number& max_prob) {
    int no_diff_obs = lambda.B.getWidth();
    int no_states = lambda.A.getWidth();

    int next_obs_guess = 0;
    number prob;
    max_prob = 0;

    //Iterate through the possible next observations
    for (int obs = 0; obs < no_diff_obs; obs++) {
        vector<number> c = vector<number>(observations.second);
        prob = 0;

        matrix* alpha_normed = hmm::a_pass(lambda, c, observations);
        for (int i = 0; i < no_states; i++) {
            number inner = 0;
            for (int j = 0; j < no_states; j++) {
                inner += lambda.A.get(j, i) * alpha_normed->get(j, observations.second - 1);
            }

            prob += lambda.B.get(i, obs) * inner;
        }

        if (prob > max_prob && !isnan(prob)) {
            next_obs_guess = obs;
            max_prob = prob;
        }

        delete alpha_normed;
    }

    return next_obs_guess;
}

//calculate the probability of seeing a certain sequence of observations from a certain HMM
number hmm::obs_seq_prob(const Lambda& lambda, const pair<vector<int>, int>& observations) {
    number log_prob = 0;

    vector<number> c = vector<number>(observations.second);

    delete hmm::a_pass(lambda, c, observations);

    for (int i = 0; i < observations.second; i++) {
        log_prob += log(c[i]);
    }

    log_prob = -log_prob;

    return log_prob;
}

//perform an alpha-pass with the given parameters and return the alpha matrix
//assumes the vector c has been initialized, and has the same length as the observation sequence
//c will be filled with the normalization constants
matrix* hmm::a_pass(const Lambda& lambda, vector<number>& c, const pair<vector<int>, int>& observations) {
    int no_states = lambda.A.getHeight();
    int seq_length = observations.second;
    
    #ifndef SAFETY_OFF_HMM
    assert(lambda.A.getWidth() == no_states);
    assert(lambda.B.getHeight() == no_states);
    assert(lambda.pi.size() == no_states);
    assert(seq_length > 0);
    assert(c.size() == seq_length);
    assert(lambda.A.row_stochastic());
    assert(lambda.B.row_stochastic());
    #endif

    c = vector<number>(seq_length);

    mat* alpha = new mat(no_states, seq_length);

    for (int i = 0; i < no_states; i++) {
        alpha->set(i, 0, lambda.pi[i] * lambda.B.get(i, observations.first[0]));
        c[0] += alpha->get(i, 0);
    }

    c[0] = 1 / c[0];

    for (int i = 0; i < no_states; i++)
        alpha->set(i, 0, alpha->get(i, 0) * c[0]);

    number sum;
    number elem;
    for (int t = 1; t < seq_length; t++) {
        for (int i = 0; i < no_states; i++) {
            sum = 0;
            for (int j = 0; j < no_states; j++)
                sum += alpha->get(j, t - 1) * lambda.A.get(j, i);
            elem = sum * lambda.B.get(i, observations.first[t]);
            alpha->set(i, t, elem);
            c[t] += elem;
        }
        c[t] = 1 / c[t];
        for (int i = 0; i < no_states; i++)
            alpha->set(i, t, alpha->get(i, t) * c[t]);
    }

    return alpha;
}

//perform a beta-pass with the given parameters and return the beta matrix
//assumes that an alpha pass has already been done so that the normalization constants in c can be used
matrix* hmm::b_pass(const Lambda& lambda, const vector<number>& c, const matrix* alpha,
                   const pair<vector<int>, int>& observations) {
    int no_states = lambda.A.getHeight();
    int seq_length = observations.second;

    #ifndef SAFETY_OFF_HMM
    assert(lambda.A.getWidth() == no_states);
    assert(lambda.B.getHeight() == no_states);
    assert(lambda.pi.size() == no_states);
    assert(seq_length > 0);
    assert(c.size() == seq_length);
    assert(alpha.getHeight() == no_states);
    assert(alpha.getWidth() == seq_length);
    assert(lambda.A.row_stochastic());
    assert(lambda.B.row_stochastic());
    #endif

    mat* beta = new mat(no_states, seq_length);

    for (int i = 0; i < no_states; i++)
        beta->set(i, seq_length - 1, c[no_states - 1]);

    for (int t = seq_length - 2; t >= 0; t--) {
        for (int i = 0; i < no_states; i++) {
            number sum = 0;
            for (int j = 0; j < no_states; j++)
                sum += lambda.A.get(i, j) * lambda.B.get(j, observations.first[t + 1]) * beta->get(j, t + 1);

            beta->set(i, t, c[t] * sum);
        }
    }

    return beta;
}

//viterbi algorithm for finding the most likely state sequence to have produced the given observation sequence
vector<int> hmm::viterbi(const Lambda& lambda, const pair<vector<int>, int>& observations) {
    int no_states = lambda.A.getHeight();
    int seq_length = observations.second;

    #ifndef SAFETY_OFF_HMM
    assert(lambda.A.getWidth() == no_states);
    assert(lambda.B.getHeight() == no_states);
    assert(lambda.pi.size() == no_states);
    assert(seq_length > 0);
    #endif

    //Use logarithms for entries to avoid underflows
    matrix log_delta = mat(no_states, seq_length);

    //To keep track of most likely sequence of states.
    vector<vector<int>> delta_index = vector<vector<int>>(no_states, vector<int>(seq_length));

    //Initialize first column
    for (int i = 0; i < no_states; i++) {
        number entry = log(lambda.B.get(i, observations.first[0])) + log(lambda.pi[i]);
        log_delta.set(i, 0, entry);
    }

    for (int t = 1; t < seq_length; t++) {
        for (int i = 0; i < no_states; i++) {

            int index_max = 0;

            //Initialize the value of max for j = 0.
            number max = log_delta.get(0, t - 1) + log(lambda.A.get(0, i));

            //Find the max_{j} of all new deltas by iterating over j.
            for (int j = 1; j < no_states; j++) {
                number new_max = log_delta.get(j, t - 1) + log(lambda.A.get(j, i));

                if (new_max > max) {
                    max = new_max;
                    index_max = j;
                }
            }

            log_delta.set(i, t, max  + log(lambda.B.get(i, observations.first[t])));
            delta_index[i][t] = index_max;
        }
    }

    vector<int> res = vector<int>(seq_length);

    //Copmute the argmax of log_deltas last column.
    number T_max = log_delta.get(0, log_delta.getWidth() - 1);
    int T_max_index;
    for (int j = 1; j < no_states; j++) {
        number new_max = log_delta.get(j, log_delta.getWidth() - 1);
        if (new_max > T_max) {
            T_max_index = j;
            T_max = new_max;
        }
    }

    res[res.size() - 1] = T_max_index;

    for (int t = res.size() - 2; t >= 0; t--) {
        res[t] = delta_index[res[t + 1]][t + 1];
    }

    return res;
}

//train the provided HMM on the given observation sequence until convergence or max_iter iterations is reached
int hmm::model_estimate(Lambda& lambda, const pair<vector<int>, int>& observations, bool verbose, int max_iter) {
    vector<number> c = vector<number>(observations.second);
    int iters = 0;
    int maxiters = max_iter;

    number old_log_prob = -std::numeric_limits<number>::infinity();

    matrix* alpha = hmm::a_pass(lambda, c, observations);
    matrix* beta = hmm::b_pass(lambda, c, alpha, observations);
    hmm::reestimate(lambda, alpha, beta, observations);

    delete alpha;
    delete beta;

    while (iters < maxiters) {
        iters++;
        number log_prob = 0;

        for (int i = 0; i < observations.second; i++) {
            log_prob += log(c[i]);
        }

        log_prob = -log_prob;

        if (verbose && ((iters & 15) == 15)) {
            cout << iters << ":\t" << log_prob << endl;
            if ((iters & 255) == 255)
                cout.flush();
        }

        if (log_prob - old_log_prob > PROB_EPSILON) {
            old_log_prob = log_prob;

            alpha = hmm::a_pass(lambda, c, observations);
            beta = hmm::b_pass(lambda, c, alpha, observations);
            hmm::reestimate(lambda, alpha, beta, observations);
        } else {
            return iters;
        }
    }

    return -iters;
}

//Alternate version of model_estimate, which uses multiple independent observation sequences
//which are thought to belong to the same model
void hmm::model_estimate(Lambda& lambda, const vector<pair<vector<int>, int>>& observations, bool verbose, int max_iter) {
    vector<number> c = vector<number>(observations[0].second);
    int iters = 0;
    int maxiters = max_iter;

    while (iters < maxiters) {
        vector<matrix*> alphas(observations.size());
        vector<matrix*> betas(observations.size());

        for (int k = 0; k < observations.size(); k++) {
            alphas[k] = hmm::a_pass(lambda, c, observations[k]);
            betas[k] = hmm::b_pass(lambda, c, alphas[k], observations[k]);
        }

        lambda = hmm::mult_seq_estimate(lambda, alphas, betas, observations);

        for (int k = 0; k < observations.size(); k++) {
            delete alphas[k];
            alphas[k] = nullptr;
            delete betas[k];
            betas[k] = nullptr;
        }

        iters++;
    }
}

//Comparison of whether two float values are within EPSILON of each other.
bool number_equal(number a, number b) {
    return abs(a - b) < EPSILON;
}

//calculates a new HMM (A, B, pi) from the data provided by an alpha and a beta pass
//A, B and pi in lambda will be overwritten
void hmm::reestimate(Lambda& lambda, const matrix* alpha, const matrix* beta, const pair<vector<int>, int>& observations) {
    int no_states = lambda.A.getHeight();
    int seq_length = observations.second;

    #ifndef SAFETY_OFF_HMM
    assert(lambda.A.getWidth() == no_states);
    assert(lambda.B.getHeight() == no_states);
    assert(lambda.pi.size() == no_states);
    assert(seq_length > 0);
    assert(alpha.getHeight() == no_states);
    assert(alpha.getWidth() == seq_length);
    assert(beta.getHeight() == no_states);
    assert(beta.getWidth() == seq_length);
    assert(lambda.A.row_stochastic());
    assert(lambda.B.row_stochastic());
    #endif

    matrix gamma = matrix(no_states, seq_length); //indexed gamma.get(i, t)
    vector<matrix*> digamma = vector<matrix*>(seq_length); //indexed digamma[t].get(i, j)
    number denom;

    for (int t = 0; t < seq_length - 1; t++) {
        denom = 0;
        digamma[t] = new matrix(no_states, no_states);

        for (int i = 0; i < no_states; i++)
            for (int j = 0; j < no_states; j++)
                denom += alpha->get(i, t) * lambda.A.get(i, j) * lambda.B.get(j, observations.first[t + 1]) * beta->get(j, t + 1);

        for (int i = 0; i < no_states; i++) {
            for (int j = 0; j < no_states; j++) {
                digamma[t]->set(i, j, alpha->get(i, t) * lambda.A.get(i, j)
                                      * lambda.B.get(j, observations.first[t + 1]) * beta->get(j, t + 1) / denom);
                gamma.set(i, t, gamma.get(i, t) + digamma[t]->get(i, j));
            }
        }
    }

    denom = 0;
    for (int i = 0; i < no_states; i++)
        denom += alpha->get(i, seq_length - 1);

    for (int i = 0; i < no_states; i++)
        gamma.set(i, seq_length - 1, alpha->get(i, seq_length - 1) / denom);

    //re-estimate A, B and pi

    //pi
    for (int i = 0; i < no_states; i++)
        lambda.pi[i] = gamma.get(i, 0);

    //A
    number numer;
    for (int i = 0; i < no_states; i++) {
        for (int j = 0; j < no_states; j++) {
            numer = 0;
            denom = 0;

            for (int t = 0; t < seq_length - 1; t++) {
                numer += digamma[t]->get(i, j);
                denom += gamma.get(i, t);
            }

            lambda.A.set(i, j, numer / denom);
        }
    }

    //B
    for (int i = 0; i < no_states; i++) {
        for (int j = 0; j < lambda.B.getWidth(); j++) {
            numer = 0;
            denom = 0;

            for (int t = 0; t < seq_length; t++) {
                if (observations.first[t] == j)
                    numer += gamma.get(i, t);
                denom += gamma.get(i, t);
            }
            lambda.B.set(i, j, numer / denom);
        }
    }

    for (int t = 0; t < seq_length; t++) {
        delete digamma[t];
        digamma[t] = nullptr;
    }
}
//calculates a new HMM (A, B, pi) from the data provided by an alpha and a beta pass
//will take into consideration multiple observation sequences
Lambda hmm::mult_seq_estimate(const Lambda& lambda, const vector<matrix*>& alphas, const vector<matrix*>& betas, const vector<pair<vector<int>, int>>& observations) {
    int no_states = lambda.A.getHeight();
    int no_observation_sequences = observations.size();

    vector<matrix*> gammas(no_observation_sequences);
    vector<vector<matrix*>*> digammas(no_observation_sequences);

    Lambda res = Lambda();

    for (int k = 0; k < no_observation_sequences; k++) {
        int seq_length = observations[k].second;

        matrix* gamma = new matrix(no_states, seq_length); //indexed gamma.get(i, t)
        vector<matrix*>* digamma = new vector<matrix*>(seq_length); //indexed digamma[t].get(i, j)
        number denom;

        for (int t = 0; t < seq_length - 1; t++) {
            denom = 0;
            (*digamma)[t] = new matrix(no_states, no_states);

            for (int i = 0; i < no_states; i++)
                for (int j = 0; j < no_states; j++)
                    denom += alphas[k]->get(i, t) * lambda.A.get(i, j) * lambda.B.get(j, observations[k].first[t + 1]) * betas[k]->get(j, t + 1);

            for (int i = 0; i < no_states; i++) {
                for (int j = 0; j < no_states; j++) {
                    (*digamma)[t]->set(i, j, alphas[k]->get(i, t) * lambda.A.get(i, j)
                                        * lambda.B.get(j, observations[k].first[t + 1]) * betas[k]->get(j, t + 1) / denom);
                    gamma->set(i, t, gamma->get(i, t) + (*digamma)[t]->get(i, j));
                }
            }
            }

        denom = 0;
        for (int i = 0; i < no_states; i++)
            denom += alphas[k]->get(i, seq_length - 1);

        for (int i = 0; i < no_states; i++)
            gamma->set(i, seq_length - 1, alphas[k]->get(i, seq_length - 1) / denom);

        digammas[k] = digamma;
        gammas[k] = gamma;
    }

    //re-estimate A, B and pi

    //pi
    for (int i = 0; i < no_states; i++) {
        res.pi[i] = 0;
        for (int k = 0; k < no_observation_sequences; k++) {
            res.pi[i] += gammas[k]->get(i, 0);
        }
        res.pi[i] /= no_observation_sequences;
    }

    //A
    number numer;
    number denom;
    for (int i = 0; i < no_states; i++) {
        for (int j = 0; j < no_states; j++) {
            numer = 0;
            denom = 0;
            for (int k = 0; k < no_observation_sequences; k++) {
                for (int t = 0; t < observations[k].second - 1; t++) {
                    numer += (*digammas[k])[t]->get(i, j);
                    denom += gammas[k]->get(i, t);
                }
            }

            res.A.set(i, j, numer / denom);
        }
    }

    //B
    for (int i = 0; i < no_states; i++) {
        for (int j = 0; j < lambda.B.getWidth(); j++) {
            numer = 0;
            denom = 0;
            for (int k = 0; k < no_observation_sequences; k++) {
                for (int t = 0; t < observations[k].second; t++) {
                    if (observations[k].first[t] == j)
                        numer += gammas[k]->get(i, t);
                    denom += gammas[k]->get(i, t);
                }
            }
            res.B.set(i, j, numer / denom);
        }
    }

    for (int k = 0; k < no_observation_sequences; k++) {
        for (int t = 0; t < observations[k].second; t++) {
            delete (*digammas[k])[t];
            (*digammas[k])[t] = nullptr;
        }
        delete digammas[k];
        digammas[k] = nullptr;
        delete gammas[k];
        gammas[k] = nullptr;
    }

    return res;
}