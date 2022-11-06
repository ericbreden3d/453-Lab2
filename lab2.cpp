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
    int pow = 2;
    int sub_n;
    int serial_result;
    double start;
    MPI_Comm cart_comm;
    MPI_Request req;
    MPI_Status stat;
    Matrix m;
    Matrix multA;
    Matrix multB;
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
        m = Matrix(n);
        m.fill_rand(1);
        m.print();
        // (m * m).print();

        if (num_procs == 1) {
            start = MPI_Wtime();
            Matrix result = m;
            for (int i = 0; i < pow; i++) {
                result = result * m;
            }
            serial_result = result.determinant();
            cout << "Serial result: " << serial_result << endl;
            cout << "Serial runtime: " << MPI_Wtime() - start << endl;
            return 0;
        }
    
        start = MPI_Wtime();
    }

    for (int i = 0; i < pow; i++) {
        if (this_rank == 0) {
            // cout << "Gettings submatrices and sending from root" << endl;
            if (i == 0) {
                multA = m;
            }
            multB = m;

            // cout << "Disassembling A" << endl;
            Matrix partsA[num_procs] = {};
            int ind = 0;
            for (int i = 0; i < n; i+=sub_n) {
                for (int j = 0; j < n; j+=sub_n) {
                    partsA[ind++] = multA.get_subm(sub_n, i, j);
                }
            }
            // cout << "Disassembling B" << endl;
            Matrix partsB[num_procs] = {};
            ind = 0;
            for (int i = 0; i < n; i+=sub_n) {
                for (int j = 0; j < n; j+=sub_n) {
                    partsB[ind++] = multB.get_subm(sub_n, i, j);
                }
            }

            // for (int i = 0; i < num_procs; i++) {
            //     partsA[i].print();
            //     partsB[i].print();
            // }

            // cout << "Distributing submatrices" << endl;
            ind = 1;
            for (int i = 0; i < dims[0]; i++) {
                for (int j = 0; j < dims[1]; j++) {
                    if (i == 0 && j == 0) continue;
                    int targ_rank;
                    int coord[2] = {i, j};
                    MPI_Cart_rank(cart_comm, coord, &targ_rank);
                    MPI_Send(partsA[ind].get_1d(), sub_n*sub_n, MPI_INT, targ_rank, 0, cart_comm);
                    MPI_Send(partsB[ind].get_1d(), sub_n*sub_n, MPI_INT, targ_rank, 0, cart_comm);
                    ind++;
                }
            }
            
            // root doesn't ned to send/recv to itself
            A = partsA[0];
            B = partsB[0];
        } else {
            // cout << "Other processes receiving" << endl;
            int buf[sub_n * sub_n];
            MPI_Recv(buf, sub_n * sub_n, MPI_INT, 0, 0, cart_comm, &stat);
            A = Matrix(buf, sub_n);
            MPI_Recv(buf, sub_n * sub_n, MPI_INT, 0, 0, cart_comm, &stat);
            B = Matrix(buf, sub_n);
        }

        // Initial Send Alignment
        // cout << "Initial alignment process" << endl;
        int A_src;
        int B_src;
        int A_dest;
        int B_dest;
        int buf[sub_n * sub_n];
        MPI_Cart_shift(cart_comm, 1, -this_coord[0], &A_src, &A_dest);
        MPI_Cart_shift(cart_comm, 0, -this_coord[1], &B_src, &B_dest);
        if (this_coord[0] != 0) {
            MPI_Send(A.get_1d(), sub_n * sub_n, MPI_INT, A_dest, 0, cart_comm);
        }
        if (this_coord[1] != 0) {
            MPI_Send(B.get_1d(), sub_n * sub_n, MPI_INT, B_dest, 0, cart_comm);
        }
        if (this_coord[0] != 0){
            MPI_Recv(buf, sub_n*sub_n, MPI_INT, A_src, 0, cart_comm, &stat);
            A = Matrix(buf, sub_n);
            cout << "built A" << endl;
        }
        if (this_coord[1] != 0) {
            MPI_Recv(buf, sub_n*sub_n, MPI_INT, B_src, 0, cart_comm, &stat);
            B = Matrix(buf, sub_n);
            cout << "built B" << endl;
        }

        MPI_Barrier(cart_comm);

        // cout << "Calculate and shift iterations" << endl;
        Matrix sum(sub_n);
        sum = sum + (A * B);
        for (int i = 1; i < dims[0]; i++) {
            MPI_Cart_shift(cart_comm, 1, -1, &A_src, &A_dest);
            MPI_Cart_shift(cart_comm, 0, -1, &B_src, &B_dest);
            MPI_Send(A.get_1d(), sub_n * sub_n, MPI_INT, A_dest, 0, cart_comm);
            MPI_Send(B.get_1d(), sub_n * sub_n, MPI_INT, B_dest, 0, cart_comm);
            MPI_Recv(buf, sub_n*sub_n, MPI_INT, A_src, 0, cart_comm, &stat);
            A = Matrix(buf, sub_n);
            MPI_Recv(buf, sub_n*sub_n, MPI_INT, B_src, 0, cart_comm, &stat);
            B = Matrix(buf, sub_n);
            sum = sum + (A * B);
        }

        // collect submatrices at root and assemble matrix
        // cout << "Collecting matrices at root" << endl;
        if (this_rank != 0) {
            MPI_Send(sum.get_1d(), sub_n * sub_n, MPI_INT, 0, 0, cart_comm);
        } else {
            Matrix parts[num_procs] = {};
            parts[0] = sum;
            int ind = 0;
            for (int i = 1; i < num_procs; i++) {
                int buf[sub_n * sub_n];
                MPI_Recv(buf, sub_n * sub_n, MPI_INT, i, 0, cart_comm, &stat);
                int coord[2];
                MPI_Cart_coords(cart_comm, i, 2, coord);
                parts[coord[1] + coord[0] * dims[0]] = Matrix(buf, sub_n);
            }
            // Matrix assem(n);
            // cout << "Assembling" << endl;
            ind = 0;
            for (int i = 0; i < n; i+=sub_n) {
                for (int j = 0; j < n; j+=sub_n) {
                    multA.add_subm(parts[ind++], sub_n, i, j);
                }
            }
            
            if (i == pow - 1) {
                cout << "Parallel result: " << multA.determinant() << endl; 
                cout << "Parallel runtime: " << MPI_Wtime() - start << endl;
                cout << endl;
            }
        }
        // assem.print();
        // cout << "Serial result: " << serial_result << endl;
        // cout << "Parallel result: " << assem.determinant() << endl;
    }

    MPI_Finalize();
}