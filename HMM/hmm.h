#ifndef AILABB_A1_hmm_H
#define AILABB_A1_hmm_H

#include "globals.h"
#include "matrix.h"

using namespace globals;

class Lambda {
    public:
        matrix A;
        matrix B;
        vector<number> pi;

        void reset();

        Lambda();
        Lambda(const matrix& A, const matrix& B, const vector<number>& pi, const vector<int> &obs_seq);
};

Lambda init_lambda();

class hmm {
    public:
        static number obs_seq_prob(const Lambda& lambda, const pair<vector<int>, int>& observations);
        static int next_obs_guess(const Lambda& lambda, const pair<vector<int>, int>& observations, number& prob);
        //static int next_obs_guess(const Lambda& lambda, pair<vector<int>, int>& observations, number& max_log_prob);
        static matrix a_pass(const Lambda& lambda, vector<number>& c, const pair<vector<int>, int>& observations);
        static matrix b_pass(const Lambda& lambda, const vector<number>& c, const matrix& alpha, const pair<vector<int>, int>& observations);
        static vector<int> viterbi(const Lambda& lambda, const pair<vector<int>, int>& observations);
        static void reestimate(Lambda& lambda, const matrix& alpha, const matrix& beta, const pair<vector<int>, int>& observations);
        static Lambda mult_seq_estimate(Lambda& lambda, const matrix& alpha, const matrix& beta, const vector<pair<vector<int>, int>>& observations);
        static void model_estimate(Lambda& lambda, const vector<pair<vector<int>, int>>& observations, bool verbose, int max_iter);
        static int model_estimate(Lambda& lambda,  const pair<vector<int>, int>& observations, bool verbose = false, int max_iter = MAX_ITERS);
};

bool number_equal(number a, number b);

#endif //AILABB_A1_hmm_H