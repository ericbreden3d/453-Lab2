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
    int sub_n;
    MPI_Comm cart_comm;
    MPI_Request req;
    MPI_Status stat;
    Matrix A;
    Matrix B;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    sub_n = n / sqrt(num_procs);

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

        Matrix parts[num_procs] = {};
        int ind = 0;
        for (int i = 0; i < n; i+=sub_n) {
            for (int j = 0; j < n; j+=sub_n) {
                // cout << i << " " << j << endl;
                parts[ind++] = m.get_subm(sub_n, j, i);
                // m.get_subm(sub_n, j, i).print();
                // parts[ind].print();
            }
        }

        ind = 1;
        for (int i = 0; i < dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {
                if (i == 0 && j == 0) continue;
                int targ_rank;
                int coord[2] = {i, j};
                MPI_Cart_rank(cart_comm, coord, &targ_rank);
                MPI_Isend(parts[ind++].get_1d(), sub_n*sub_n, MPI_INT, targ_rank, 0, cart_comm, &req);
            }
        }
        
        // root doesn't ned to send/recv to itself
        A = parts[0];
        B = parts[0];
    } else {
        int buf[sub_n * sub_n];
        MPI_Recv(buf, sub_n * sub_n, MPI_INT, 0, 0, cart_comm, &stat);
        A = Matrix(buf, sub_n);
        B = A;
    }

    // Initial Alignment
    for (int i = 0; i < dims[0]; i++) {
        for (int j = 0; j < dims[1]; j++) {
            if (this_coord[0] == i && this_coord[1] == j) {
                int A_src;
                int B_src;
                int A_dest;
                int B_dest;
                if (i != 0) {
                    MPI_Cart_shift(cart_comm, 0, i, &A_src, &A_dest);
                    MPI_Isend(A.get_1d(), sub_n * sub_n, MPI_INT, A_dest, 0, cart_comm, &req);
                    int buf[sub_n*sub_n];
                    MPI_Recv(buf, sub_n*sub_n, MPI_INT, A_src, 0, cart_comm, &stat);
                    A = Matrix(buf, sub_n);

                    int A_coord[2];
                    MPI_Cart_coords(cart_comm, A_dest, 2, A_coord);
                    // cout << i << "," << j << " " << "sending to " << A_coord[0] << "," << A_coord[1] << endl;

                    cout << "DONE A " << this_coord[0] << "," << this_coord[1] << endl;

                }
                if (j != 0) {
                    MPI_Cart_shift(cart_comm, 1, j, &B_src, &B_dest);
                    MPI_Isend(B.get_1d(), sub_n * sub_n, MPI_INT, B_dest, 0, cart_comm, &req);
                    int buf[sub_n*sub_n];
                    MPI_Recv(buf, sub_n*sub_n, MPI_INT, B_src, 0, cart_comm, &stat);
                    B = Matrix(buf, sub_n);
                    cout << "DONE B " << this_coord[0] << "," << this_coord[1] << endl;
                    int B_coord[2];
                    MPI_Cart_coords(cart_comm, B_dest, 2, B_coord);
                    // cout << i << "," << j << " " << "sending to " << B_coord[0] << "," << B_coord[1] << endl;
                }
                
            }
        }
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
