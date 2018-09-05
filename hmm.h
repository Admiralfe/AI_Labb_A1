#ifndef AILABB_A1_hmm_H
#define AILABB_A1_hmm_H

#include "globals.h"
#include "matrix.h"

using namespace globals;

class hmm {
    public:
        static number a_pass(const matrix& A, const matrix& B, const vector<number>& pi, const vector<number>& obs_seq);
};

#endif //AILABB_A1_hmm_H