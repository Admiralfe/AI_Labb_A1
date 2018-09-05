#include <iostream>
#include <cmath>
#include <vector>
#include "matrix.h"
#include "globals.h"
#include "hmm.h"

typedef vector<number> vec;
typedef matrix mat;

int main() {
    /*
    mat A(3, 2);
    A.fill({
        1, 2,
        3, 5,
        7, 9
    });

    vec x = {
        1,
        2
    };

    mat B(2, 4);

    B.fill({
        1, 2, 3, 4,
        5, 6, 7, 8
    });



    vec y = A * x;
    mat C = A * B;

    cout << A << "times\n" << B << "is\n" << C << endl;

    cout << y << endl;

    cout << A.get_col(1);

    return 0;
     */

    mat transition = matrix::parse_stdin();
    mat emission = matrix::parse_stdin();
    mat state = matrix::parse_stdin();
    vec state_vec = vec(state.getWidth());
    for (int i = 0; i < state.getWidth(); i++)
        state_vec[i] = state.get(0, i);

    vec observation_sequence = matrix::parse_vec_stdin();
    vec normalization = vec(observation_sequence.size());

    vec alphat = hmm::a_pass(transition, emission, state_vec, observation_sequence, normalization);

    number sum = 0;
    for (int i = 0; i < alphat.size(); i++)
        sum += alphat[i];
    sum /= normalization[normalization.size() - 1];

    cout << sum << endl;
    return 0;
}