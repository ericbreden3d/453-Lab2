#include <iostream>
#include <string>
#include <math.h>
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
    int n = stoi(argv[1]);
    int sub_sz;
    MPI_Comm cart_comm;
    MPI_Request req;
    MPI_Status stat;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    sub_sz = n / sqrt(num_procs);

    MPI_Dims_create(num_procs, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &cart_comm);

    // get new rank, cart coords, and amount of procs in each dim
    MPI_Comm_rank(cart_comm, &this_rank);
    MPI_Cart_coords(cart_comm, this_rank, 2, this_coord);
    get_dim_counts(2, cart_comm, dim_counts);

    if (this_rank == 0) {
        Matrix m(n);
        m.fill_rand();
        m.print();

        int* parts[num_procs];
        int ind = 0;
        for (int i = 0; i < n; i+=sub_sz) {
            for (int j = 0; j < n; j+=sub_sz) {
                // cout << i << " " << j << endl;
                Matrix sub = m.get_subm(sub_sz, i, j);
                parts[ind] = sub.get_1d();
                for (int k = 0; k < sub_sz * sub_sz; k++) {
                    cout << parts[ind][k] << " ";
                }
                ind++;
            }
            cout << endl;
        }

        for (int i = 0; i < num_procs; i++) {
            for (int j = 0; j < sub_sz * sub_sz; j++) {
                cout << parts[i][j] << " ";
            }
            cout << endl;
        }

        // MPI_Isend(m.get_1d(), n*n, MPI_INT, 1, 0,
        //       cart_comm, &req);
        // cout << "Determinant: " << m.determinant() << endl;
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
