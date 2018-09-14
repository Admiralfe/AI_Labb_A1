#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

#include "hmm.h"
#include "globals.h"
#include "matrix.h"

#define NO_OBS 9
#define NO_HS 6
#define TIME_OUT -2
#define MAX_ITERS_REACHED -1

using namespace globals;
using namespace std;

typedef matrix mat;
typedef vector<number> vec;

Lambda::Lambda() {
    A = matrix::random_uniform(NO_HS, NO_HS, 0.1);
    B = matrix::random_uniform(NO_HS, NO_OBS, 0.1);
    pi = matrix::random_uniform(1, NO_HS, 0.1).get_row(0);

    obs_seq = vector<int>(100);
    no_obs = 0;
}

Lambda::Lambda(const matrix& transition, const matrix& emission, const vector<number>& init_state, const vector<int>& observations) {
    A = transition;
    B = emission;
    pi = init_state;
    obs_seq = observations;
}

//gör en alpha-pass med givna parameterar och returnerar alpha-matrisen,
//förutsätter att vektorn c är initialiserad med nollor och har samma längd som obs_seq
//c kommer att populeras med normeringskonstanter
matrix hmm::a_pass(const Lambda& lambda, vector<number>& c) {
    int no_states = lambda.A.getHeight();
    int seq_length = lambda.no_obs;

    //cerr << "XD " << seq_length << endl;
    
    assert(lambda.A.getWidth() == no_states);
    assert(lambda.B.getHeight() == no_states);
    assert(lambda.pi.size() == no_states);
    assert(seq_length > 0);
    assert(c.size() == seq_length);
    assert(lambda.A.row_stochastic());
    assert(lambda.B.row_stochastic());

    c = vector<number>(seq_length);

    mat alpha = mat(no_states, seq_length);

    for (int i = 0; i < no_states; i++) {
        alpha.set(i, 0, lambda.pi[i] * lambda.B.get(i, lambda.obs_seq[0]));
        c[0] += alpha.get(i, 0);
    }

    c[0] = 1 / c[0];

    for (int i = 0; i < no_states; i++)
        alpha.set(i, 0, alpha.get(i, 0) * c[0]);

    number sum;
    number elem;
    for (int t = 1; t < seq_length; t++) {
        for (int i = 0; i < no_states; i++) {
            sum = 0;
            for (int j = 0; j < no_states; j++)
                sum += alpha.get(j, t - 1) * lambda.A.get(j, i);
            elem = sum * lambda.B.get(i, lambda.obs_seq[t]);
            alpha.set(i, t, elem);
            c[t] += elem;
        }
        c[t] = 1 / c[t];
        for (int i = 0; i < no_states; i++)
            alpha.set(i, t, alpha.get(i, t) * c[t]);
    }

    return alpha;
}

//förutsätter att en alpha-pass redan gjorts och att värdena i c inte har ändrats
matrix hmm::b_pass(const Lambda& lambda, const vector<number>& c, const matrix& alpha) {
    int no_states = lambda.A.getHeight();
    int seq_length = lambda.no_obs;

    assert(lambda.A.getWidth() == no_states);
    assert(lambda.B.getHeight() == no_states);
    assert(lambda.pi.size() == no_states);
    assert(seq_length > 0);
    assert(c.size() == seq_length);
    assert(alpha.getHeight() == no_states);
    assert(alpha.getWidth() == seq_length);
    assert(lambda.A.row_stochastic());
    assert(lambda.B.row_stochastic());

    mat beta = mat(no_states, seq_length);

    for (int i = 0; i < no_states; i++)
        beta.set(i, seq_length - 1, c[no_states - 1]);

    for (int t = seq_length - 2; t >= 0; t--) {
        for (int i = 0; i < no_states; i++) {
            number sum = 0;
            for (int j = 0; j < no_states; j++)
                sum += lambda.A.get(i, j) * lambda.B.get(j, lambda.obs_seq[t + 1]) * beta.get(j, t + 1);

            beta.set(i, t, c[t] * sum);
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

vector<int> hmm::viterbi(const Lambda& lambda) {
    int no_states = lambda.A.getHeight();
    int seq_length = lambda.no_obs;

    assert(lambda.A.getWidth() == no_states);
    assert(lambda.B.getHeight() == no_states);
    assert(lambda.pi.size() == no_states);
    assert(seq_length > 0);

    //Use logarithms for entries to avoid underflows
    matrix log_delta = mat(no_states, seq_length);

    //To keep track of most likely sequence of states.
    vector<vector<int>> delta_index = vector<vector<int>>(no_states, vector<int>(seq_length));

    //Initialize first column
    for (int i = 0; i < no_states; i++) {
        number entry = log(lambda.B.get(i, lambda.obs_seq[0])) + log(lambda.pi[i]);
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

            log_delta.set(i, t, max  + log(lambda.B.get(i, lambda.obs_seq[t])));
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

int hmm::model_estimate(Lambda& lambda, const Deadline& pDue, bool verbose, int max_iter) {
    vector<number> c = vector<number>(lambda.no_obs);
    int iters = 0;
    int maxiters = max_iter;

    number old_log_prob = -std::numeric_limits<double>::infinity();

    cerr << lambda.A << endl;

    matrix alpha = hmm::a_pass(lambda, c);
    matrix beta = hmm::b_pass(lambda, c, alpha);
    hmm::reestimate(lambda, alpha, beta);

    while (iters < maxiters) {
        if (pDue.remainingMs() < 5) {
            return TIME_OUT;
        }
        iters++;
        number log_prob = 0;

        for (int i = 0; i < lambda.no_obs; i++) {
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

            alpha = hmm::a_pass(lambda, c);
            beta = hmm::b_pass(lambda, c, alpha);
            hmm::reestimate(lambda, alpha, beta);
        } else {
            return iters;
        }
    }

    return MAX_ITERS_REACHED;
}

//Comparison if two float values are within EPSILON of each other.
bool number_equal(number a, number b) {
    return abs(a - b) < EPSILON;
}

//räknar ut ny modell (A, B, pi) utifrån datan från en alpha- och en beta-pass
//A, B och pi kommer att skrivas över
void hmm::reestimate(Lambda& lambda, const matrix& alpha, const matrix& beta) {
    int no_states = lambda.A.getHeight();
    int seq_length = lambda.no_obs;

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

    matrix gamma = matrix(no_states, seq_length); //indexed gamma.get(i, t)
    vector<matrix*> digamma = vector<matrix*>(seq_length); //indexed digamma[t].get(i, j)
    number denom;

    for (int t = 0; t < seq_length - 1; t++) {
        denom = 0;
        digamma[t] = new matrix(no_states, no_states);

        for (int i = 0; i < no_states; i++)
            for (int j = 0; j < no_states; j++)
                denom += alpha.get(i, t) * lambda.A.get(i, j) * lambda.B.get(j, lambda.obs_seq[t + 1]) * beta.get(j, t + 1);

        for (int i = 0; i < no_states; i++) {
            for (int j = 0; j < no_states; j++) {
                digamma[t]->set(i, j, alpha.get(i, t) * lambda.A.get(i, j) * lambda.B.get(j, lambda.obs_seq[t + 1]) * beta.get(j, t + 1) / denom);
                gamma.set(i, t, gamma.get(i, t) + digamma[t]->get(i, j));
            }
        }
    }

    denom = 0;
    for (int i = 0; i < no_states; i++)
        denom += alpha.get(i, seq_length - 1);

    for (int i = 0; i < no_states; i++)
        gamma.set(i, seq_length - 1, alpha.get(i, seq_length - 1) / denom);

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
                if (lambda.obs_seq[t] == j)
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