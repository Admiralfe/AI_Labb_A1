#ifndef AILABB_A1_hmm_H
#define AILABB_A1_hmm_H

#include "globals.h"
#include "matrix.h"

using namespace globals;

typedef struct {
    vector<number> pi;
    matrix transition;
    matrix emission;
} HM_model;

class hmm {
    public:
        static matrix a_pass(const matrix& A, const matrix& B, const vector<number>& pi, const vector<int>& obs_seq, vector<number>& c);
        static matrix b_pass(const matrix& A, const matrix& B, const vector<number>& pi, const vector<int>& obs_seq, const vector<number>& c, const matrix& alpha);
        static vector<int> viterbi(matrix A, matrix B, vector<number> pi, vector<int> obs_seq);
        static void reestimate(matrix& A, matrix& B, vector<number> pi, const vector<int>& obs_seq, const matrix& alpha, const matrix& beta);
        HM_model model_estimate(matrix A, matrix B, vector<number> pi) {

        }
};

#endif //AILABB_A1_hmm_H