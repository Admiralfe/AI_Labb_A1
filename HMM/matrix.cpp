#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <set>

#include "matrix.h"
#include "globals.h"
#include "hmm.h"

using namespace globals;

//#define SAFETY_OFF_MATRIX

matrix::matrix(int height, int width){
    #ifdef SAFETY_OFF_MATRIX
    assert(height >= 1 && width >= 1);
    #endif
    this->height = height;
    this->width = width;
    this->elements = vector<vector<number>>(height, vector<number>(width));
    
    /*
    for (int i = 0; i < height; i++)
        this->elements[i] = vector<T>(width);
        */
}

vector<number> matrix::operator *(const vector<number>& v) {
    #ifdef SAFETY_OFF_MATRIX
    assert(v.size() == this->width);
    #endif

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
    #ifdef SAFETY_OFF_MATRIX
    assert(this->width == m.height);
    #endif

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

matrix matrix::operator *(number x) {
    matrix res(height, width);
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
            res.elements[i][j] = this->elements[i][j] * x;
    
    return res;
}
matrix matrix::operator +(const matrix& m) {
    #ifdef SAFETY_OFF_MATRIX
    assert(m.height == this->height);
    assert(m.width == this->width);
    #endif

    matrix res(height, width);
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
            res.elements[i][j] = this->elements[i][j] + m.elements[i][j];
    
    return res;
}
matrix matrix::operator -(const matrix& m) {
    #ifdef SAFETY_OFF_MATRIX
    assert(m.height == this->height);
    assert(m.width == this->width);
    #endif

    matrix res(height, width);
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
            res.elements[i][j] = this->elements[i][j] - m.elements[i][j];
    
    return res;
}


bool matrix::row_stochastic() const {
    for (int i = 0; i < this->height; i++) {
        number row_sum = 0;

        for (int j = 0; j < this->width; j++) {
            if (this->get(i, j) < 0)
                return false;
            row_sum += this->get(i, j);
        }

        if (!number_equal(row_sum, 1))
            return false;
    }

    return true;
}

number matrix::get(int i, int j) const {
    #ifdef SAFETY_OFF_MATRIX
    assert(i >= 0 && i < height && j >= 0 && j < width);
    #endif
    return this->elements[i][j];
}

void matrix::set(int i, int j, number element) {
    #ifdef SAFETY_OFF_MATRIX
    assert(i >= 0 && i < height && j >= 0 && j < width);
    #endif
    this->elements[i][j] = element;
}

void matrix::fill(const vector<number>& v) {
    #ifdef SAFETY_OFF_MATRIX
    assert(v.size() == width * height);
    #endif

    for (int i = 0; i < height; i++)
        for (int j = 0, offset = i * width; j < width; j++)
            this->elements[i][j] = v[offset + j];
}

matrix matrix::random_uniform(int h, int w, number variance) {
    matrix res = matrix(h, w);

    long seed = std::chrono::system_clock::now().time_since_epoch().count();

    default_random_engine generator(seed);
    uniform_real_distribution<number> distribution(-variance, variance);

    number uniform_value = 1 / (number) w;

    for (int i = 0; i < h; i++) {
        number current_row_sum = 0;
        for (int j = 0; j < w; j++) {
            //We set the last element in each row so that the matrix becomes row stochastic,
            //that is each row must sum to 1.
            if (j == (w - 1)) {
                //Make sure the last element doesn't become negative.
                while (current_row_sum > 1) {
                    //Divide all previous elements by 1.1 until we no longer have negative last value.
                    for (int k = 0; k < w -1; k++) {
                        res.set(i, k, res.get(i, k) / 1.1);
                    }

                    current_row_sum /= 1.1;
                }

                res.set(i, j, 1 - current_row_sum);
            }
            else {
                //We don't want negative values in a row stochastic matrix
                number entry = max(0.0, uniform_value + distribution(generator));
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
    res << fixed << setprecision(FLOAT_PRECISION);

    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            res << elements[row][col] << "\t";
        }
        res << endl;
    }

    return res.str();
}

vector<number> matrix::get_col(int j) const {
    vector<number> res = vector<number>(height);

    for (int i = 0; i < height; i++)
        res[i] = elements[i][j];

    return res;
}

vector<number> matrix::get_row(int i) const {
    vector<number> res = vector<number>(width);

    for (int j = 0; j < width; j++)
        res[j] = elements[i][j];

    return res;
}

void matrix::to_stdout() const {
    cout << this->height << ' ' << this->width;

    for (int row = 0; row < this->height; row++) {
        for (int col = 0; col < this->width; col++) {
            cout << ' ' << this->get(row, col);
        }
    }
}

//use -1 for infinity norm
number matrix::distance(const matrix& other, int norm) const {
    #ifdef SAFETY_OFF_MATRIX
    assert(this->height == other.height);
    assert(this->width == other.width);
    #endif

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

number matrix::row_distance_squared(const matrix& other, int i1, int i2) const {
    number sum = 0;
    for (int j = 0; j < width; j++)
        sum += (elements[i1][j] - other.elements[i2][j]) * (elements[i1][j] - other.elements[i2][j]);

    return sum;
}

//beräkna det ungefärliga kvadrerade avståndet mellan två matriser, utan hänsyn till radnumrering
number matrix::distance_squared(const matrix& other) const {
    #ifdef SAFETY_OFF_MATRIX
    assert(height == other.height);
    assert(width == other.width);
    #endif

    number total_dist = 0;

    std::set<int> left_to_check;
    for (int i = 0; i < height; i++)
        left_to_check.insert(i);

    for (int i = 0; i < height; i++) {
        number mindist = numeric_limits<number>::max();
        int minrow = -1;
        number dist = mindist;

        for (auto row : left_to_check) {
            dist = this->row_distance_squared(other, i, row);
            if (dist < mindist) {
                mindist = dist;
                minrow = row;
            }
        }

        total_dist += mindist;
        left_to_check.erase(minrow);
    }

    return total_dist;
}

//beräkna det ungefärliga kvadrerade avståndet mellan två matriser, utan hänsyn till radnumrering
//sparar radnumreringen för översättning mellan matrisernas index i reordering, vars .size() ska vara lika med this->height
//om square_reordered_matrices så används istället elementen i reordering för att jämföra de två matriserna, som ska vara kvadratiska
number matrix::distance_squared(const matrix& other, vector<int>& reordering, bool square_reordered_matrices) const {
    #ifdef SAFETY_OFF_MATRIX
    assert(height == other.height);
    assert(width == other.width);
    assert(reordering.size() == height);
    #endif

    number total_dist = 0;

    if (square_reordered_matrices) {
        #ifdef SAFETY_OFF_MATRIX
        assert(width == height);
        for (auto index : reordering)
            assert(index >= 0 && index < width);
        #endif
        
        number x;
        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++) {
                x = (elements[i][j] - other.elements[reordering[i]][reordering[j]]);
                total_dist += x * x;
            }
    } else {
        std::set<int> left_to_check;
        for (int i = 0; i < height; i++)
            left_to_check.insert(i);

        for (int i = 0; i < height; i++) {
            number mindist = numeric_limits<number>::max();
            int minrow = -1;
            number dist = mindist;

            for (auto row : left_to_check) {
                dist = this->row_distance_squared(other, i, row);
                if (dist < mindist) {
                    mindist = dist;
                    minrow = row;
                }
            }

            reordering[i] = minrow;

            total_dist += mindist;
            left_to_check.erase(minrow);
        }
    }

    return total_dist;
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
    #ifdef SAFETY_OFF_MATRIX
    assert(a.size() == b.size());
    #endif

    number sum = 0;
    for (int i = 0; i < a.size(); i++)
        sum += a[i] * b[i];
    
    return sum;
}