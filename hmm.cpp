#include <cassert>
#include <iostream>
#include <cmath>

#include "hmm.h"
#include "globals.h"
#include "matrix.h"

using namespace globals;
using namespace std;

typedef matrix mat;
typedef vector<number> vec;

vector<number> hmm::a_pass(const matrix& A, const matrix& B, const vector<number>& pi, const vector<int>& obs_seq, vector<number>& c) {
    int no_states = A.getHeight();
    int seq_length = obs_seq.size();

    assert(A.getWidth() == no_states);
    assert(B.getHeight() == no_states);
    assert(pi.size() == no_states);
    assert(seq_length > 0);
    assert(c.size() == seq_length);

    mat alpha = mat(no_states, seq_length);

    for (int i = 0; i < no_states; i++) {
        alpha.set(i, 0, pi[i] * B.get(i, obs_seq[0]));
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
                sum += alpha.get(j, t - 1) * A.get(j, i);
            elem = sum * B.get(i, obs_seq[t]);
            alpha.set(i, t, elem);
            c[t] += elem;
        }
        c[t] = 1 / c[t];
        for (int i = 0; i < no_states; i++)
            alpha.set(i, t, alpha.get(i, t) * c[t]);
    }

    return alpha.get_col(seq_length - 1);
}

vector<int> hmm::viterbi(matrix A, matrix B, vector<number> pi, vector<int> obs_seq) {
    int no_states = A.getHeight();
    int seq_length = obs_seq.size();

    assert(A.getWidth() == no_states);
    assert(B.getHeight() == no_states);
    assert(pi.size() == no_states);
    assert(seq_length > 0);

    //Use logarithms for entries to avoid underflows
    matrix log_delta = mat(no_states, seq_length);

    //To keep track of most likely sequence of states.
    vector<vector<int>> delta_index = vector<vector<int>>(no_states, vector<int>(seq_length));

    //Initialize first column
    for (int i = 0; i < no_states; i++) {
        number entry = log(B.get(i, obs_seq[0])) + log(pi[i]);
        log_delta.set(i, 0, entry);
    }


    int index_max;
    //Loops over states
    for (int i = 0; i < no_states; i++) {
        //Loops over time steps
        for (int t = 1; t < seq_length; t++) {

            //Initialize the value of max for j = 0.
            number max = log_delta.get(t - 1, 0) + log(A.get(0, i)) + log(B.get(i, obs_seq[0]));

            //Find the max_{j} of all new deltas by iterating over j.
            for (int j = 1; j < no_states; j++) {
                number new_max = log_delta.get(t - 1, j) + log(A.get(j, i)) + log(B.get(i, obs_seq[j]));

                if (new_max > max) {
                    max = new_max;
                    index_max = j;
                }
            }

            log_delta.set(i, t, max);
            delta_index[i][t] = index_max;
        }
    }

    vector<int> res = vector<int>(seq_length);

    //index_max should still have the relevant max value of the last time step since the for loop finished.
    int max_state_index = index_max;
    res[res.size() - 1] = max_state_index;

    for (int t = res.size() - 2; t >= 0; t--) {
        res[t] = delta_index[max_state_index][t];
        max_state_index = delta_index[max_state_index][t];
    }

    return res;
}