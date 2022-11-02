#include <iostream>
#include <stdlib.h>
#include <time.h>
using namespace std;

void make_matrix(int n) {
    srand(time(NULL));
    int** m = new int*[n];
    for (int i = 0; i < n; i++) {
        cout << m[i] << endl;
    }
}

int main(int argc, int** argv) {
    int this_rank;
    int num_procs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (this_rank == 0) {
        make_matrix(10);
    }
}


// free