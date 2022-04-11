#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[]) {
    int rank, num_proc; 

    MPI_Init(&argc, &argv); // initialize the MPI environment
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc); // get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get the rank of the process

    cout << "Hello world from rank " << rank << endl;

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}