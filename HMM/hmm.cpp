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
#define NO_HS 4
#define TIME_OUT -2
//#define ALWAYS_ROW_STOCHASTIC
#define SAFETY_OFF_HMM

using namespace globals;
using namespace std;

typedef matrix mat;
typedef vector<number> vec;

Lambda::Lambda() {
    A = matrix::random_uniform(NO_HS, NO_HS, 0.1);
    B = matrix::random_uniform(NO_HS, NO_OBS, 0.1);
    pi = matrix::random_uniform(1, NO_HS, 0.1).get_row(0);
}

Lambda::Lambda(const matrix& transition, const matrix& emission, const vector<number>& init_state, const vector<int>& observations) {
    A = transition;
    B = emission;
    pi = init_state;
}

//Resets the model parameters of lambda to random_uniform values but keeps the observation sequence.
void Lambda::reset() {
    A = matrix::random_uniform(NO_HS, NO_HS, 0.1);
    B = matrix::random_uniform(NO_HS, NO_OBS, 0.1);
    pi = matrix::random_uniform(1, NO_HS, 0.1).get_row(0);
}

/*
 * Returns the most likely next observation given a model lambda and a current observation sequence.
 * max_log_prob will be set to the probability of that most likely observation upon function return.
 */
/*int hmm::next_obs_guess(const Lambda& lambda, pair<vector<int>, int>& observations, number& max_log_prob) {
    int no_diff_obs = lambda.B.getWidth();

    number norm_factor = 0;
    max_log_prob = -std::numeric_limits<number>::infinity();
    number log_prob;
    int next_obs_guess = 0;

    //cerr << "In next_obs_guess: " << lambda.obs_seq[lambda.no_obs - 1] << endl;

    //Iterate through the possible next observations
    for (int obs = 0; obs < no_diff_obs; obs++) {
        vector<number> c = vector<number>(observations.second + 1);

        //Add next observation guess to the sequence
        observations.first[observations.second] = obs;
        observations.second++;
        matrix alpha_normed = hmm::a_pass(lambda, c, observations);
        //Compute log probability and pick maximum probability one.
        log_prob = 0;
        for (int i = 0; i < observations.second; i++) {
            log_prob += log(c[i]);
        }
        log_prob = -log_prob;

        //cerr << "current (log) probability: " << log_prob << endl;
        //cerr << "max (log) probability: " << max_log_prob << endl;
        if (log_prob > max_log_prob) {
            next_obs_guess = obs;
            max_log_prob = log_prob;
        }

        norm_factor += exp(log_prob);

        //This effectively removes the added observation guess from obs_seq.
        observations.second--;
    }

    return next_obs_guess;
}
 */

int hmm::next_obs_guess(const Lambda& lambda, const pair<vector<int>, int>& observations, number& max_prob) {
    int no_diff_obs = lambda.B.getWidth();
    int no_states = lambda.A.getWidth();

    int next_obs_guess = 0;
    number prob;
    max_prob = 0;

    //cerr << "In next_obs_guess: " << lambda.obs_seq[lambda.no_obs - 1] << endl;

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

            //cerr << "inner sum: " << inner << endl;
            prob += lambda.B.get(i, obs) * inner;
        }

        if (prob > max_prob && !isnan(prob)) {
            next_obs_guess = obs;
            max_prob = prob;
        }

        delete alpha_normed;

        //This effectively removes the added observation guess from obs_seq.
        //lambda.no_obs--;
    }

    return next_obs_guess;
}

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
//gör en alpha-pass med givna parameterar och returnerar alpha-matrisen,
//förutsätter att vektorn c är initialiserad med nollor och har samma längd som obs_seq
//c kommer att populeras med normeringskonstanter
matrix* hmm::a_pass(const Lambda& lambda, vector<number>& c, const pair<vector<int>, int>& observations) {
    int no_states = lambda.A.getHeight();
    int seq_length = observations.second;

    //cerr << "XD " << seq_length << endl;
    
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

//förutsätter att en alpha-pass redan gjorts och att värdena i c inte har ändrats
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

    /*
    int current_max;
    int current;

    vector<int> res = vector<int>(seq_length);

    for (int t = 0; t < seq_length; t++) {
        current_max = -1;
        for (int i = 0; i < no_states; i++) {
            current = alpha.get(i, t) * beta.get(i, t);
            if (current > current_max) {
                current_max = current;
                res[t] = i;
            }
        }
    }*/

    return beta;
}

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

int hmm::model_estimate(Lambda& lambda, const pair<vector<int>, int>& observations, bool verbose, int max_iter) {
    vector<number> c = vector<number>(observations.second);
    int iters = 0;
    int maxiters = max_iter;

    number old_log_prob = -std::numeric_limits<number>::infinity();

    matrix* alpha = hmm::a_pass(lambda, c, observations);
    matrix* beta = hmm::b_pass(lambda, c, alpha, observations);
    hmm::reestimate(lambda, alpha, beta, observations);


    while (iters < maxiters) {
        iters++;
        number log_prob = 0;

        for (int i = 0; i < observations.second; i++) {
            log_prob += log(c[i]);
        }

        log_prob = -log_prob;

        if (verbose && ((iters & 15) == 15)) {
            //cout << iters << ":\t" << log_prob << endl;
            //if ((iters & 255) == 255)
            //    cout.flush();
        }

        if (log_prob - old_log_prob > PROB_EPSILON) {
            old_log_prob = log_prob;

            delete alpha;
            delete beta;

            alpha = hmm::a_pass(lambda, c, observations);
            beta = hmm::b_pass(lambda, c, alpha, observations);
            hmm::reestimate(lambda, alpha, beta, observations);

            //cerr << lambda.A << endl;
        } else {
            return iters;
        }
    }

    return -iters;
}

/*
 * Overloaded version of model_estimate, which uses multiple independent observation sequences which are thought to
 * belong to the same model.
 */
void hmm::model_estimate(Lambda& lambda, const vector<pair<vector<int>, int>>& observations, bool verbose, int max_iter) {
    vector<number> c = vector<number>(observations[0].second);
    int iters = 0;
    int maxiters = max_iter;
    int no_states = lambda.A.getWidth();
    int no_obs = lambda.B.getWidth();

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

//Comparison if two float values are within EPSILON of each other.
bool number_equal(number a, number b) {
    return abs(a - b) < EPSILON;
}

//räknar ut ny modell (A, B, pi) utifrån datan från en alpha- och en beta-pass
//A, B och pi kommer att skrivas över
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

    #ifdef ALWAYS_ROW_STOCHASTIC
    bool A_needed = false, B_needed = false, pi_needed = false;
    //A
    for (int i = 0; i < no_states; i++) {
        number sum = 0;
        for (int j = 0; j < no_states; j++)
            sum += lambda.A.get(i, j);
        if (!number_equal(sum, 1)) {
            A_needed = true;
            cerr << "A_" << i << " needed fixing (" << sum << ")" << endl;
            sum = 1 / sum;
            for (int j = 0; j < no_states; j++)
                lambda.A.set(i, j, lambda.A.get(i, j) * sum);
        }
    }
    //B
    for (int i = 0; i < no_states; i++) {
        number sum = 0;
        for (int j = 0; j < lambda.B.getWidth(); j++)
            sum += lambda.B.get(i, j);
        if (!number_equal(sum, 1)) {
            B_needed = true;
            cerr << "B_" << i << " needed fixing (" << sum << ")" << endl;
            sum = 1 / sum;
            for (int j = 0; j < lambda.B.getWidth(); j++)
                lambda.B.set(i, j, lambda.B.get(i, j) * sum);
        }
    }
    //pi
    {
        number sum = 0;
        for (int j = 0; j < lambda.pi.size(); j++)
            sum += lambda.pi[j];
        if (!number_equal(sum, 1)) {
            pi_needed = true;
            cerr << "pi needed fixing (" << sum << ")" << endl;
            sum = 1 / sum;
            for (int j = 0; j < lambda.pi.size(); j++)
                lambda.pi[j] *= sum;
        }
    }

    if (A_needed)
        cerr << "A is now" << endl << lambda.A << endl;
        
    if (B_needed)
        cerr << "B is now" << endl << lambda.B << endl;
        
    if (pi_needed)
        cerr << "pi is now" << endl << lambda.pi << endl;
    #endif
}

Lambda hmm::mult_seq_estimate(const Lambda& lambda, const vector<matrix*>& alphas, const vector<matrix*>& betas, const vector<pair<vector<int>, int>>& observations) {
    int no_states = lambda.A.getHeight();
    int no_observation_sequences = observations.size();

    vector<matrix*> gammas(no_observation_sequences);
    vector<vector<matrix*>*> digammas(no_observation_sequences);

    Lambda res = Lambda();

    for (int k = 0; k < no_observation_sequences; k++) {
        int seq_length = observations[k].second;

        matrix gamma = matrix(no_states, seq_length); //indexed gamma.get(i, t)
        vector<matrix*> digamma = vector<matrix*>(seq_length); //indexed digamma[t].get(i, j)
        number denom;

        for (int t = 0; t < seq_length - 1; t++) {
            denom = 0;
            digamma[t] = new matrix(no_states, no_states);

            for (int i = 0; i < no_states; i++)
                for (int j = 0; j < no_states; j++)
                    denom += alphas[k]->get(i, t) * lambda.A.get(i, j) * lambda.B.get(j, observations[k].first[t + 1]) * betas[k]->get(j, t + 1);

            for (int i = 0; i < no_states; i++) {
                for (int j = 0; j < no_states; j++) {
                    digamma[t]->set(i, j, alphas[k]->get(i, t) * lambda.A.get(i, j)
                                        * lambda.B.get(j, observations[k].first[t + 1]) * betas[k]->get(j, t + 1) / denom);
                    gamma.set(i, t, gamma.get(i, t) + digamma[t]->get(i, j));
                }
            }
        }

        denom = 0;
        for (int i = 0; i < no_states; i++)
            denom += alphas[k]->get(i, seq_length - 1);

        for (int i = 0; i < no_states; i++)
            gamma.set(i, seq_length - 1, alphas[k]->get(i, seq_length - 1) / denom);
        
        digammas[k] = &digamma;
        gammas[k] = &gamma;
    }

    //re-estimate A, B and pi

    //pi
    for (int i = 0; i < no_states; i++) {
        res.pi[i] = 0;
        for (int k = 0; k < no_observation_sequences; k++)
            res.pi[i] += gammas[k]->get(i, 0);
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
#ifdef ALWAYS_ROW_STOCHASTIC
    bool A_needed = false, B_needed = false, pi_needed = false;
    //A
    for (int i = 0; i < no_states; i++) {
        number sum = 0;
        for (int j = 0; j < no_states; j++)
            sum += lambda.A.get(i, j);
        if (!number_equal(sum, 1)) {
            A_needed = true;
            cerr << "A_" << i << " needed fixing (" << sum << ")" << endl;
            sum = 1 / sum;
            for (int j = 0; j < no_states; j++)
                lambda.A.set(i, j, lambda.A.get(i, j) * sum);
        }
    }
    //B
    for (int i = 0; i < no_states; i++) {
        number sum = 0;
        for (int j = 0; j < lambda.B.getWidth(); j++)
            sum += lambda.B.get(i, j);
        if (!number_equal(sum, 1)) {
            B_needed = true;
            cerr << "B_" << i << " needed fixing (" << sum << ")" << endl;
            sum = 1 / sum;
            for (int j = 0; j < lambda.B.getWidth(); j++)
                lambda.B.set(i, j, lambda.B.get(i, j) * sum);
        }
    }
    //pi
    {
        number sum = 0;
        for (int j = 0; j < lambda.pi.size(); j++)
            sum += lambda.pi[j];
        if (!number_equal(sum, 1)) {
            pi_needed = true;
            cerr << "pi needed fixing (" << sum << ")" << endl;
            sum = 1 / sum;
            for (int j = 0; j < lambda.pi.size(); j++)
                lambda.pi[j] *= sum;
        }
    }

    if (A_needed)
        cerr << "A is now" << endl << lambda.A << endl;

    if (B_needed)
        cerr << "B is now" << endl << lambda.B << endl;

    if (pi_needed)
        cerr << "pi is now" << endl << lambda.pi << endl;
#endif
}