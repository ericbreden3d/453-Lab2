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
        // Matrix a(n);
        // Matrix b(n);
        // a.fill_rand();
        // b.fill_rand();
        // Matrix res = a + b;
        // a.print();
        // b.print();
        // res.print();
        // cout << "dete: " << res.determinant();
        // return 0;

        Matrix m(n);
        m.fill_rand();
        m.print();

        Matrix parts[num_procs] = {};
        int ind = 0;
        for (int i = 0; i < n; i+=sub_n) {
            for (int j = 0; j < n; j+=sub_n) {
                // cout << i << " " << j << endl;
                parts[ind++] = m.get_subm(sub_n, i, j);
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

    // if (this_rank == 1) {
        // cout << this_coord[0] << "," << this_coord[1] << endl;
        // A.print();
    // }
        

    // Initial Send Alignment
    int A_src;
    int B_src;
    int A_dest;
    int B_dest;
    int buf[sub_n * sub_n];
    MPI_Cart_shift(cart_comm, 1, this_coord[0], &A_src, &A_dest);
    MPI_Cart_shift(cart_comm, 0, this_coord[1], &B_src, &B_dest);
    if (this_coord[0] != 0) {
        MPI_Isend(A.get_1d(), sub_n * sub_n, MPI_INT, A_dest, 0, cart_comm, &req);
    }
    if (this_coord[1] != 0) {
        MPI_Isend(B.get_1d(), sub_n * sub_n, MPI_INT, B_dest, 0, cart_comm, &req);
    }
    if (this_coord[0] != 0){
        MPI_Recv(buf, sub_n*sub_n, MPI_INT, A_src, 0, cart_comm, &stat);
        A = Matrix(buf, sub_n);
    }
    if (this_coord[1] != 0) {
        MPI_Recv(buf, sub_n*sub_n, MPI_INT, B_src, 0, cart_comm, &stat);
        B = Matrix(buf, sub_n);
    }

    // if (this_rank == 2) {
    //     A.print();
    //     B.print();
    // }

    MPI_Barrier(cart_comm);
    cout << "Initial alignment complete\n";

    // Calc and shift
    Matrix sum(sub_n);
    sum = sum + (A * B);
    // if (this_rank == 1) {
    //     A.print();
    //     B.print();
    //     sum.print();
    // }
    for (int i = 1; i < dims[0]; i++) {
        MPI_Cart_shift(cart_comm, 1, 1, &A_src, &A_dest);
        MPI_Cart_shift(cart_comm, 0, 1, &B_src, &B_dest);
        MPI_Isend(A.get_1d(), sub_n * sub_n, MPI_INT, A_dest, 0, cart_comm, &req);
        MPI_Isend(B.get_1d(), sub_n * sub_n, MPI_INT, B_dest, 0, cart_comm, &req);
        MPI_Recv(buf, sub_n*sub_n, MPI_INT, A_src, 0, cart_comm, &stat);
        A = Matrix(buf, sub_n);
        MPI_Recv(buf, sub_n*sub_n, MPI_INT, B_src, 0, cart_comm, &stat);
        B = Matrix(buf, sub_n);
        sum = sum + (A * B);
    }

    // if (this_rank == 1)
    //     sum.print();
    
    if (this_rank == 0) {
        Matrix parts[num_procs] = {};
        parts[0] = sum;
        for (int i = 1; i < num_procs; i++) {
            int buf[sub_n * sub_n];
            MPI_Recv(buf, sub_n * sub_n, MPI_INT, i, 0, cart_comm, &stat);
            cout << i << "received" << endl;
            int coord[2];
            MPI_Cart_coords(cart_comm, i, 2, coord);
            parts[coord[0] + coord[1] * dims[0]] = Matrix(buf, sub_n);
            parts[i].print();
        }
    } else {
        MPI_Isend(sum.get_1d(), sub_n * sub_n, MPI_INT, 0, 0, cart_comm, &req);
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
