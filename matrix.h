

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
    matrix operator *(number x);
    matrix operator +(const matrix& m);
    matrix operator -(const matrix& m);
    
    matrix(int height, int width);
    number get(int i, int j) const;
    void set(int i, int j, number element);
    void fill(const vector<number>& v);

    int getHeight() const { return height; }
    int getWidth() const { return width; }

    string to_string() const;
    vector<number> get_col(int j) const;
    vector<number> get_row(int i) const;

    static matrix parse_stdin(); 
    static vector<int> parse_intvec_stdin();
    static matrix random_uniform(int h, int w, double variance);

    void to_stdout() const;

    bool row_stochastic() const;
    number distance(const matrix& other, int norm = 2) const;

private:
    int width, height;
    vector<vector<number>> elements;
};

ostream& operator<< (ostream& outs, const matrix& m);

ostream& operator<< (ostream& outs, const vector<number>& v);

ostream& operator<< (ostream& outs, const vector<int>& v);

number operator *(const vector<number>& a, const vector<number>& b);

#endif //AILABB_A1_matrix_H
