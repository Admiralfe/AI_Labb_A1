#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <chrono>

#include "matrix.h"
#include "globals.h"
#include "hmm.h"

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

    number sum;

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

bool matrix::row_stochastic() const {
    for (int i = 0; i < this->height; i++) {
        number row_sum = 0;

        for (int j = 0; j < this->width; j++) {
            row_sum += this->get(i, j);
        }

        if (!number_equal(row_sum, 1))
            return false;
    }

    return true;
}

number matrix::get(int i, int j) const {
    assert(i >= 0 && i < height && j >= 0 && j < width);
    return this->elements[i][j];
}

void matrix::set(int i, int j, number element) {
    assert(i >= 0 && i < height && j >= 0 && j < width);
    this->elements[i][j] = element;
}

void matrix::fill(const vector<number>& v) {
    assert(v.size() == width * height);

    for (int i = 0; i < height; i++)
        for (int j = 0, offset = i * width; j < width; j++)
            this->elements[i][j] = v[offset + j];
}

matrix matrix::random_uniform(int h, int w, double variance) {
    matrix res = matrix(h, w);

    long seed = std::chrono::system_clock::now().time_since_epoch().count();

    default_random_engine generator(seed);
    uniform_real_distribution<double> distribution(-variance, variance);

    double uniform_value = 1 / (double) w;

    for (int i = 0; i < h; i++) {
        double current_row_sum = 0;
        for (int j = 0; j < w; j++) {
            //We set the last element in each row so that the matrix becomes row stochastic,
            //that is each row must sum to 1.
            if (j == (w - 1))
                res.set(i, j, 1 - current_row_sum);
            else {
                double entry = uniform_value + distribution(generator);
                current_row_sum += entry;
                res.set(i, j, entry);
            }
        }
    }

    return res;
}

//Parses a matrix given on standard input in the format:
//[number_of_rows number_of_columns entries] separated by spaces.
// Each line on standard input corresponds to a new matrix.
matrix matrix::parse_stdin() {
    //Will contain stdin split by whitespace.
    vector<string> input;
    vector<number> number_vec;

    string line;
    getline(cin, line);

    stringstream stream;
    stream.str(line);

    string s;

    while (getline(stream, s, ' ')) {
        input.push_back(s);
    }


    int height = stoi(input[0]);
    int width = stoi(input[1]);

    matrix res = matrix(height, width);

    for (int i = 2; i < input.size(); i++) {
        number_vec.push_back(stod(input[i]));
    }

    res.fill(number_vec);

    return res;
}

//se ovan
vector<int> matrix::parse_intvec_stdin() {
    //Will contain stdin split by whitespace.
    vector<string> input;

    string line;
    getline(cin, line);

    stringstream stream;
    stream.str(line);

    string s;

    while (getline(stream, s, ' ')) {
        input.push_back(s);
    }

    int length = stoi(input[0]);

    vector<int> res = vector<int>(length);

    for (int i = 0; i < length; i++)
        res[i] = stoi(input[i+1]);

    return res;
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

void matrix::to_stdout() {
    cout << this->height << ' ' << this->width;

    for (int row = 0; row < this->height; row++) {
        for (int col = 0; col < this->width; col++) {
            cout << ' ' << this->get(row, col);
        }
    }
}

//use -1 for infinity norm
number matrix::distance(const matrix& other, int norm) const {
    assert(this->height == other.height);
    assert(this->width == other.width);

    number res = 0;
    switch(norm) {
        case 1:
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    res += abs(elements[i][j] - other.elements[i][j]);
            break;
        
        case 2:
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    res += pow(elements[i][j] - other.elements[i][j], 2);

            res = sqrt(res);
            break;

        case -1: //infinity
            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    res = max(res, abs(elements[i][j] - other.elements[i][j]));
            break;

        default:
            throw new invalid_argument("norm");
    }

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

ostream& operator<< (ostream& outs, const vector<int>& v) {
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