#ifndef AILABB_A1_hmm_H
#define AILABB_A1_hmm_H

#include "globals.h"
#include "matrix.h"

using namespace globals;

class hmm {
    public:
        static matrix a_pass(const matrix& A, const matrix& B, const vector<number>& pi, const vector<int>& obs_seq, vector<number>& c);
        static matrix b_pass(const matrix& A, const matrix& B, const vector<number>& pi, const vector<int>& obs_seq, const vector<number>& c, const matrix& alpha);
        static vector<int> viterbi(matrix A, matrix B, vector<number> pi, vector<int> obs_seq);
};

#endif //AILABB_A1_hmm_H