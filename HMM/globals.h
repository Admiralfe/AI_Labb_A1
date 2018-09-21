
#ifndef AILABB_A1_GLOBALS
#define AILABB_A1_GLOBALS

namespace globals {
    //define the type used in most calculations. can be changed to float or long double, for instance
    typedef double number;
    //floating point comparison threshold
    const number EPSILON = 0.00001;
    //convergence limit
    const number PROB_EPSILON = 0.00001;
    //maximum number of iterations of model training unless (usually) overriden in the function call
    const int MAX_ITERS = 30;
    //number of decimals to use when printing matrices
    const int FLOAT_PRECISION = 3;
}

#endif //AILABB_A1_GLOBALS
