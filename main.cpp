#include <iostream>
#include <cmath>
#include <vector>
#include "matrix.h"
#include "globals.h"

typedef vector<number> vec;
typedef matrix mat;

int main() {
    /*
    mat A(3, 2);
    A.fill({
        1, 2,
        3, 5,
        7, 9
    });

    vec x = {
        1,
        2
    };

    mat B(2, 4);

    B.fill({
        1, 2, 3, 4,
        5, 6, 7, 8
    });

    vec y = A * x;
    mat C = A * B;

    cout << A << "times\n" << B << "is\n" << C << endl;

    for (int i = 0; i < y.size(); i++)
        std::cout << y[i] << std::endl;
    return 0;
     */

    matrix m = matrix::parse_stdin();

    cout << m;
}