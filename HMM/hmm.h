#ifndef AILABB_A1_hmm_H
#define AILABB_A1_hmm_H

#include "globals.h"
#include "matrix.h"
#include "../duckhunt/Deadline.hpp"

using namespace globals;
using namespace ducks;

class Lambda {
    public :
        matrix A;
        matrix B;
        vector<number> pi;
        vector<int> obs_seq;
        int no_obs;

        Lambda();
        Lambda(const matrix& A, const matrix& B, const vector<number>& pi, const vector<int> &obs_seq);
        
        static vector<int> group_models(const vector<Lambda>& hmms, int number_of_groups, bool best_guess = false);
};

Lambda init_lambda();

class hmm {
    public:
        static matrix a_pass(const Lambda& lambda, vector<number>& c);
        static matrix b_pass(const Lambda& lambda, const vector<number>& c, const matrix& alpha);
        static vector<int> viterbi(const Lambda& lambda);
        static void reestimate(Lambda& lambda, const matrix& alpha, const matrix& beta);
        static int model_estimate(Lambda& lambda,  const Deadline& pDue, bool verbose = false, int max_iter = MAX_ITERS);
};

bool number_equal(number a, number b);

#endif //AILABB_A1_hmm_H