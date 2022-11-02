#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <mpi.h>
using namespace std;

class Matrix {
    int** matrix;
    int size;

 public:
    Matrix(int n) {
        size = n;
        int** m = new int*[size];
        for (int i = 0; i < size; i++) {
            m[i] = new int[size];
            for (int j = 0; j < size; j++) {
                m[j][i] = 3 - rand() % 4;
            }
        }
    }

    void print_matrix() {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                cout << matrix[j][i] << " ";
            }
            cout << endl;
        }
    }
};

int main(int argc, char** argv) {
    int this_rank;
    int num_procs;
    int n = stoi(argv[1])

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (this_rank == 0) {
        m = Matrix(n);
        m.print_matrix();
    }
}


// free