

#ifndef AILABB_A1_matrix_H
#define AILABB_A1_matrix_H

#include <vector>
#include <string>
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
    vector<number> get_col(int j);

    static matrix parse_stdin();

    void to_stdout();

private:
    int width, height;
    vector<vector<number>> elements;
};

ostream& operator<< (ostream& outs, const matrix& m);

ostream& operator<< (ostream& outs, const vector<number>& v);

number operator *(const vector<number>& a, const vector<number>& b);

#endif //AILABB_A1_matrix_H
