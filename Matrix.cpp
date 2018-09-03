#include <cassert>
#include <iostream>

#include "Matrix.h"
#include "globals.h"

using namespace globals;

vector<number> Matrix::operator *(vector<number>& v) {
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

Matrix::Matrix(int height, int width){
    this->height = height;
    this->width = width;
    this->elements = vector<vector<number>>(height, vector<number>(width));
    
    /*
    for (int i = 0; i < height; i++)
        this->elements[i] = vector<T>(width);
        */
}

number Matrix::get(int i, int j) {
    return this->elements[i][j];
}

void Matrix::set(int i, int j, number element) {
    this->elements[i][j] = element;
}