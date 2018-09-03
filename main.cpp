#include <iostream>
#include <cmath>
#include <vector>
#include "Matrix.h"


int main() {
    std::cout << "Hello, World!" << std::endl;
    Matrix<double> A(3, 2);
    vector<double> x(2);

    A.set(0, 0, 1);
    A.set(0, 1, 2);
    A.set(1, 0, 5);
    A.set(1, 1, 7);
    A.set(2, 0, 13);
    A.set(2, 1, 17);

    x[0] = 1;
    x[1] = 2;


    vector<double> y = A * x;
    std::cout << y.size() << std::endl;



    for (int i = 0; i < y.size(); i++)
        std::cout << y[i] << std::endl;
    return 0;
}