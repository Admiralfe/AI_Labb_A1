#ifndef AILABB_A1_hmm_H
#define AILABB_A1_hmm_H

#include "globals.h"
#include "matrix.h"

using namespace globals;

typedef struct Lambda {
    matrix A;
    matrix B;
    vector<number> pi;
    vector<int> obs_seq;
};

Lambda init_lambda();

class hmm {
    public:
        static matrix a_pass(const matrix& A, const matrix& B, const vector<number>& pi, const vector<int>& obs_seq, vector<number>& c);
        static matrix b_pass(const matrix& A, const matrix& B, const vector<number>& pi, const vector<int>& obs_seq, const vector<number>& c, const matrix& alpha);
        static vector<int> viterbi(matrix A, matrix B, vector<number> pi, vector<int> obs_seq);
        static void reestimate(matrix& A, matrix& B, vector<number>& pi, const vector<int>& obs_seq, const matrix& alpha, const matrix& beta);
        static int model_estimate(matrix& A, matrix& B, vector<number>& pi, vector<int> obs_seq, bool verbose = false, int max_iter = MAX_ITERS);
};

bool number_equal(number a, number b);

vector<int> lambda_group_models(const vector<Lambda>& hmms, int number_of_groups);

#endif //AILABB_A1_hmm_H