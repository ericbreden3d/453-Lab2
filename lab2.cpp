#include <iostream>
#include <string>
// #include <mpi.h>
#include "Matrix.h"
using namespace std;


int main(int argc, char** argv) {
    int this_rank;
    int num_procs;
    int n = stoi(argv[1]);

    // MPI_Init(&argc, &argv);
    // MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    this_rank = 0;
    if (this_rank == 0) {
        
        
    }

    Matrix m(n);
    Matrix j(n);
    m.fill_rand();
    // for (int i = 0; i < INT_MAX/2; i++) {}  // delay
    j.fill_rand();
    // m.print();
    Matrix w = m*j;
    // m.print();
    // j.print();
    w.print();

    // Matrix sub = w.get_subm(4, 3, 2);
    w.determinant();
}



// MPI_Wtime()
