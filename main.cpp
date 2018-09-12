#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

#include "matrix.h"
#include "globals.h"
#include "hmm.h"

typedef vector<number> vec;
typedef matrix mat;

void mainC() {
    vector<int> obs_seq = matrix::parse_intvec_stdin();
    
    mat A = mat(3, 3);
    A.fill({
        0.54, 0.26, 0.20,
        0.19, 0.53, 0.28,
        0.22, 0.18, 0.6
    });

    mat B = mat(3, 4);
    B.fill({
        0.5, 0.2, 0.11, 0.19,
        0.22, 0.28, 0.23, 0.27,
        0.19, 0.21, 0.15, 0.45
    });

    vec pi = {0.3, 0.2, 0.5};

    int iterations = hmm::model_estimate(A, B, pi, obs_seq, true, 5000);

    cout << "A" << endl << A << endl;
    cout << "B" << endl << B << endl;
    cout << "pi" << endl << pi << endl;

    if (iterations != -1)
        cout << "Converged after " << iterations << " iterations" << endl;
    else
        cout << "Did not converge" << endl;
}

void mainDE() {
    mat transition = matrix::parse_stdin();
    mat emission = matrix::parse_stdin();
    mat state = matrix::parse_stdin();
    vec state_vec = vec(state.getWidth());
    
    for (int i = 0; i < state.getWidth(); i++)
        state_vec[i] = state.get(0, i);

    vector<int> observation_sequence = matrix::parse_intvec_stdin();

    assert(transition.row_stochastic());
    assert(emission.row_stochastic());

    hmm::model_estimate(transition, emission, state_vec, observation_sequence);

    transition.to_stdout();
    cout << endl;
    emission.to_stdout();
    cout << endl;
}

void mainTest() {
    matrix matrix = matrix::random_uniform(5, 10, 0.01);

    cout << matrix << endl;
}

int main() {
    //mainDE();
    //mainC();
    mainTest();
    return 0;
}