#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <math.h>
#include "Matrix.h"
using namespace std;

Matrix Matrix::get_detrm_subm(int split) {
    // cout << "Getting sub-matrix" << endl;
    Matrix new_m(size - 1);
    for (int i = 0; i < this->size; i++) {
        if (i == split) {
            split = -1;
            continue;
        }
        int x = split == -1 ? i - 1 : i;  // if already split, shift left
        for (int j = 1; j < this->size; j++) {
            new_m(x, j - 1) = (*this)(i, j);
            arr[x + (j - 1) * size] = (*this)(i, j);
        }
    }
    return new_m;
}

int Matrix::detrm_helper(Matrix m) {
    if (m.size == 1) {
        return m(0, 0);
    }
    int op = 1;
    int detrm = 0;
    for (int i = 0; i < m.size; i++) {
        int rec = detrm_helper(m.get_detrm_subm(i));
        int result = m(i, 0) * rec;
        if (op) {
            detrm += result;
            op--;
        } else {
            detrm -= result;
            op++;
        }
    }
    return detrm;
}

Matrix::Matrix(int n) {
    // cout << "Creating matrix" << endl;
    size = n;
    matrix = new int*[size];
    arr = new int[size * size];
    for (int i = 0; i < size; i++) {
        matrix[i] = new int[size];
        for (int j = 0; j < size; j++) {
            matrix[i][j] = 0;
            arr[i + j * size] = 0;
        }
    }
}

Matrix::Matrix(int* buff, int len) {
    size = sqrt(len);
    matrix = new int*[size];
    for (int i = 0; i < size; i++) {
        matrix[i] = new int[size];
    }
    arr = new int[len];
    int i = 0;
    int j = 0;
    for (int ind = 0; ind < len; ind++) {
        if (ind == size) {
            j++;
            i = 0;
        }
        matrix[i][j] = buff[ind];
        arr[i + j * size] = buff[ind];
    }
}

Matrix::~Matrix() {
    // cout << "Deleting matrix" << endl;
    for (int i = 0; i < size; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
    delete[] arr;
}

void Matrix::fill_rand() {
    // cout << "Filling with random values" << endl;
    srand(time(NULL));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matrix[i][j] = 3 - rand() % 4;
            arr[i + j * size] = 0;
        }
    }
}

Matrix Matrix::operator*(Matrix& other) {
    Matrix new_m(size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int sum = 0;
            for (int k = 0; k < size; k++) {
                // cout << (*this)(k, i) << " " << other(j, k) << endl;
                // cout << (*this)(i, k) << " " << other(k, j) << endl;
                sum += (*this)(k, i) * other(j, k);
            }
            new_m(j, i) = sum;
        }
    }
    return new_m;
}

Matrix Matrix::get_subm(int len, int x, int y) {
    Matrix new_m(len);
    for (int i = x; i < len + x; i++) {
        for (int j = y; j < len + y; j++) {
            new_m(i - x, j - y) = (*this)(i, j);
            new_m.arr[i - x + (j - y) * len] = (*this)(i, j);
        }
    }
    return new_m;
}

int* Matrix::get_1d() {
    return arr;
}

int Matrix::determinant() {
    int result = detrm_helper(*this);
    cout << "Determinant : " << result << endl;
    return result;
}

int& Matrix::operator()(int x, int y) {
    return matrix[x][y];
}

void Matrix::print() {
    cout << "Printing matrix" << endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            cout << (*this)(j, i) << " ";
        }
        cout << endl;;
    }
    cout << endl;
}