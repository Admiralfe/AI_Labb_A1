#include <cassert>
#include <iostream>
#include <sstream>

#include "matrix.h"
#include "globals.h"

using namespace globals;

matrix::matrix(int height, int width){
    assert(height >= 1 && width >= 1);
    this->height = height;
    this->width = width;
    this->elements = vector<vector<number>>(height, vector<number>(width));
    
    /*
    for (int i = 0; i < height; i++)
        this->elements[i] = vector<T>(width);
        */
}

vector<number> matrix::operator *(const vector<number>& v) {
    assert(v.size() == this->width);

    vector<number> res(this->height);

    int sum;

    for (int row = 0; row < this->height; row++) {
        sum = 0;
        for (int col = 0; col < this->width; col++) {
            sum += this->elements[row][col] * v[col];
        }
        res[row] = sum;
    }

    return res;
}

matrix matrix::operator*(const matrix& m) {
    assert(this->width == m.height);

    matrix res(this->height, m.width);

    int sum;

    for (int col = 0; col < m.width; col++) {
        for (int row = 0; row < this->height; row++) {
            sum = 0;
            for (int term = 0; term < this->width; term++) {
                sum += this->elements[row][term] * m.elements[term][col];
            }

            res.elements[row][col] = sum;
        }
    }

    return res;
}

number matrix::get(int i, int j) {
    assert(0 >= i && i < height && j >= 0 && j < width);
    return this->elements[i][j];
}

void matrix::set(int i, int j, number element) {
    assert(0 >= i && i < height && j >= 0 && j < width);
    this->elements[i][j] = element;
}

void matrix::fill(const vector<number>& v) {
    assert(v.size() == width * height);

    for (int i = 0; i < height; i++)
        for (int j = 0, offset = i * width; j < width; j++)
            elements[i][j] = v[offset + j];
}

string matrix::to_string() const {
    stringstream res;

    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            res << elements[row][col] << "\t";
        }
        res << endl;
    }

    return res.str();
}

vector<number> matrix::get_col(int j) {
    vector<number> res = vector<number>(height);

    for (int i = 0; i < height; i++)
        res[i] = elements[i][j];

    return res;
}

ostream& operator<< (ostream& outs, const matrix& m) {
	return outs << m.to_string();
}

ostream& operator<< (ostream& outs, const vector<number>& v) {
    for (int i = 0; i < v.size(); i++)
	    outs << v[i] << endl;
    return outs;
}

number operator *(const vector<number>& a, const vector<number>& b) {
    assert(a.size() == b.size());

    number sum = 0;
    for (int i = 0; i < a.size(); i++)
        sum += a[i] * b[i];
    
    return sum;
}