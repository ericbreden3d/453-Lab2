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
    int neighbors[4] = {};
    MPI_Comm cart_comm;
    int n = stoi(argv[1]);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    MPI_Dims_create(num_procs, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &cart_comm);

    // get new rank, cart coords, and amount of procs in each dim
    MPI_Comm_rank(cart_comm, &this_rank);
    MPI_Cart_coords(cart_comm, this_rank, 2, this_coord);
    get_dim_counts(2, cart_comm, dim_counts);

    if (this_rank == 0) {
        Matrix m(n);
        m.fill_rand();
        MPI_Request req;
        MPI_Isend(m.get_1d(), n*n, MPI_INT, 1, 0,
              cart_comm, &req);
    }

    if (this_rank == 1) {
        MPI_Status stat;
        int buf[n];
        MPI_Recv(buf, n, MPI_INT, 0, 0,
             cart_comm, stat);
        for (int i = 0; i < n; i++) {
            cout << buf[i] << " ";
        }
        cout << endl;
    }
}




    
    // if (this_rank == 0) {
    //     int nl_coord[2] = {3, 0};
    //     int nl_rank;
    //     int nr_coord[2] = {1, 0};
    //     int nr_rank;
    //     int nu_coord[2] = {0, 1};
    //     int nu_rank;
    //     int nd_coord[2] = {0, 3};
    //     int nd_rank;
    //     MPI_Cart_shift(cart_comm, 0, 1, &neighbors[0], &neighbors[1]);
    //     MPI_Cart_shift(cart_comm, 1, 1, &neighbors[2], &neighbors[3]);
    //     MPI_Cart_rank(cart_comm, nl_coord, &nl_rank);
    //     MPI_Cart_rank(cart_comm, nr_coord, &nr_rank);
    //     MPI_Cart_rank(cart_comm, nu_coord, &nu_rank);
    //     MPI_Cart_rank(cart_comm, nd_coord, &nd_rank);
    //     cout << "Dims: " << dim_counts[0] << " " << dim_counts[1] << endl;
    //     cout << "This coord: " << this_coord[0] << "," << this_coord[1] << endl;
    //     cout << "This rank: " << this_rank << endl;
    //     cout << "Left neighbor rank " << nl_rank << endl;
    //     cout << "Right neighbor rank " << nr_rank << endl;
    //     cout << "Up neighbor rank " << nu_rank << endl;
    //     cout << "Down neighbor rank " << nd_rank << endl;
    //     cout << "Neighbors:";
    //     for (int n : neighbors) {
    //         cout << " " << n;
    //     }
    //     cout << endl;  



// MPI_Wtime()
