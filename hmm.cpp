#include <cassert>
#include "hmm.h"
#include "globals.h"
#include "matrix.h"

using namespace globals;

typedef matrix mat;
typedef vector<number> vec;

vector<number> hmm::a_pass(const matrix& A, const matrix& B, const vector<number>& pi, const vector<number>& obs_seq) {
    int no_states = A.getHeight();
    int seq_length = obs_seq.size();
    
    assert(A.getWidth() == no_states);
    assert(B.getHeight() == no_states);
    assert(pi.size() == no_states);
    assert(seq_length > 0);

    mat alpha = mat(no_states, seq_length);
    vec c = vec(seq_length);

    for (int i = 0; i < no_states; i++) {
        alpha.set(i, 0, pi[i] * B.get(i, obs_seq[0]));
        c[0] += alpha.get(i, 0);
    }

    c[0] = 1 / c[0];

    for (int i = 0; i < no_states; i++)
        alpha.set(i, 0, alpha.get(i, 0) * c[0]);

    int sum;
    int elem;
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