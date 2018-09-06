#include <iostream>
#include <cmath>
#include <vector>
#include "matrix.h"
#include "globals.h"
#include "hmm.h"

typedef vector<number> vec;
typedef matrix mat;

int main() {
    cout << "oh look it compiled";
    cout.flush();


    mat transition = matrix::parse_stdin();
    mat emission = matrix::parse_stdin();
    mat state = matrix::parse_stdin();
    vec state_vec = vec(state.getWidth());

    for (int i = 0; i < state.getWidth(); i++)
        state_vec[i] = state.get(0, i);

    vector<int> observation_sequence = matrix::parse_intvec_stdin();

    vec normalization = vec(observation_sequence.size());

    //mat alpha = hmm::a_pass(transition, emission, state_vec, observation_sequence, normalization);;

    //vector<int> most_likely_states = hmm::b_pass(transition, emission, state_vec, observation_sequence, normalization, alpha);

    vector<int> most_likely_states = hmm::viterbi(transition, emission, state_vec, observation_sequence);
    for (int i = 0; i < most_likely_states.size(); i++)
        cout << most_likely_states[i] << " ";
        
    cout << endl;
    return 0;
}