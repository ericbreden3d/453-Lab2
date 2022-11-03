#include <iostream>
#include <string>
#include <mpi.h>
#include "Matrix.h"
using namespace std;

// get # of processors in each dim - dims_arr is out param
void get_dim_counts(int m, MPI_Comm comm, int* dims_arr) {
    int dim_num;
    int periods[m];
    int coords[m];
    MPI_Cart_get(comm, m, dims_arr, periods, coords);
}

int main(int argc, char** argv) {
    int this_rank;
    int num_procs;
    int dims[2] = {0, 0};
    int dim_counts[2];
    int periods[2] = {true, true};
    int this_coord[2];
    int neighbors[4];
    MPI_Comm = cart_comm;
    int n = stoi(argv[1]);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    MPI_Dims_create(size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &cart_comm);

    // get new rank, cart coords, and amount of procs in each dim
    MPI_Comm_rank(new_comm, &this_rank);
    MPI_Cart_coords(new_comm, this_rank, m, this_coord);
    get_dim_counts(m, new_comm, dim_counts);


    if (this_rank == 0) {
        MPI_Cart_shift(cart_comm, 0, 1, &neighbours[LEFT], &neighbours);
        cout <"This rank: " << this_rank << endl;
        cout << "Neighbors:";
        for (int n : neighbors) {
            cout << " n";
        }
        cout << endl;
    }
}



// MPI_Wtime()
