

#ifndef AILABB_A1_MATRIX_H
#define AILABB_A1_MATRIX_H

#include <vector>
#include "globals.h"

using namespace std;
using namespace globals;

class Matrix {
public:
    vector<number> operator *(vector<number>& v);
    
    Matrix(int height, int width);
    number get(int i, int j);
    void set(int i, int j, number element);

private:
    int width, height;
    vector<vector<number>> elements;
};


#endif //AILABB_A1_MATRIX_H
