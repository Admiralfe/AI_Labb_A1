#ifndef AILABB_A1_classification_H
#define AILABB_A1_classification_H

#include <vector>
#include "globals.h"
#include "hmm.h"

class classification {
    public:
        static std::vector<int> group_models(const std::vector<Lambda>& hmms, int number_of_groups, bool best_guess = false);
};

#endif //AILABB_A1_classification_H