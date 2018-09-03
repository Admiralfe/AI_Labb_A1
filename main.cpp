#include <iostream>
#include <cmath>
#include <vector>
#include "Matrix.h"
#include "globals.h"

typedef vector<number> vec;
typedef Matrix mat;

int main() {
    std::cout << "test\n";
    mat A(3, 2);
    vec x(2);

    cout << "initialized\n";
    cout.flush();

    A.set(0, 0, 1);
    A.set(0, 1, 2);
    A.set(1, 0, 5);
    A.set(1, 1, 7);
    A.set(2, 0, 13);
    A.set(2, 1, 17);

    cout << "assigned A\n";


    x[0] = 1;
    x[1] = 2;

    vec y = A * x;
    std::cout << y.size() << std::endl;

    for (int i = 0; i < y.size(); i++)
        std::cout << y[i] << std::endl;
    return 0;
}