

#ifndef AILABB_A1_MATRIX_H
#define AILABB_A1_MATRIX_H

#include <vector>
#include <cassert>
#include <iostream>

using namespace std;

template<class T>
class Matrix {
public:
    vector<T> operator *(vector<T>& v) {
        assert(v.size() == this->width);

        vector<T> res(this->height);

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
    Matrix(int height, int width) {
        this->height = height;
        this->width = width;
        this->elements = vector<vector<T>>(height, vector<T>(width));
        /*
        for (int i = 0; i < height; i++)
            this->elements[i] = vector<T>(width);
            */
    }
    T get(int i, int j) {
        return this->elements[i][j];
    }
    void set(int i, int j, T element) {
        this->elements[i][j] = element;
    }

private:
    int width, height;
    vector<vector<T>> elements;
};


#endif //AILABB_A1_MATRIX_H
