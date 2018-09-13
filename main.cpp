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
    /*
    mat A = mat(3, 3);
    A.fill({
        0.69, 0.04, 00.27,
        0.08, 0.77, 0.15,
        0.18, 0.34, 0.48
    });

    assert(A.row_stochastic());

    mat B = mat(3, 4);
    B.fill({
        0.67, 0.25, 0.06, 0.02,
        0.13, 0.37, 0.26, 0.24,
        0.05, 0.13, 0.15, 0.67
    });

    vec pi = {1, 0, 0};


    mat Actual = mat(3, 3);
    Actual.fill({
                   0.7, 0.05, 0.25,
                   0.1, 0.8, 0.1,
                   0.2, 0.3, 0.5
           });

    mat Bctual = mat(3, 4);
    Bctual.fill({
                   0.7, 0.2, 0.1, 0,
                   0.1, 0.4, 0.3, 0.2,
                   0, 0.1, 0.2, 0.7
           });

    cout << A << endl;
    cout << B << endl;
    */

    mat A = matrix::random_uniform(2, 2, 0.05);
    mat B = matrix::random_uniform(2, 4, 0.05);

    cout << A << endl;
    cout << B << endl;

    vec pi = matrix::random_uniform(1, 2, 0.05).get_row(0);



    int iterations = hmm::model_estimate(A, B, pi, obs_seq, true, 5000);

    cout << "A" << endl << A << endl;
    cout << "B" << endl << B << endl;
    cout << "pi" << endl << pi << endl;

    if (iterations != -1)
        cout << "Converged after " << iterations << " iterations" << endl;
    else
        cout << "Did not converge" << endl;

    //cout << A - Actual << endl;
    //cout << B - Bctual << endl;


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
    matrix matrix = matrix::random_uniform(3, 3, 0.01);

    cout << matrix << endl;
}

int old_main() {
    //mainDE();
    mainC();
    //mainTest();
    return 0;
}