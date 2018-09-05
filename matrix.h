

#ifndef AILABB_A1_matrix_H
#define AILABB_A1_matrix_H

#include <vector>
#include "globals.h"

using namespace std;
using namespace globals;

class matrix {
public:
    vector<number> operator *(const vector<number>& v);
    matrix operator *(const matrix& m);
    
    matrix(int height, int width);
    number get(int i, int j);
    void set(int i, int j, number element);
    void fill(const vector<number>& v);

    string to_string() const;

    static matrix parse_stdin();

private:
    int width, height;
    vector<vector<number>> elements;
};

ostream& operator<< (ostream& outs, const matrix& m);


#endif //AILABB_A1_matrix_H
