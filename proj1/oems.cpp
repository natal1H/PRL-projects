#include <iostream>
#include <fstream>
#include <mpi.h>

#define ROOT 0
#define LEN_NUMBERS 8

void print_processor(int rank, int num_proc, int low, int high) {
    // print processor info
    std::cout << "Rank " << rank << "/" << num_proc << ": low=" << low << ", high=" << high << std::endl;
}

int main(int argc, char *argv[]) {
    int neighbour_matrix[19][2] = {
        {4, 5},   // rank 0
        {4, 5},   // rank 1
        {6, 7},   // rank 2
        {6, 7},   // rank 3
        {10, 8},  // rank 4
        {8, 13},  // rank 5
        {10, 9},  // rank 6
        {9, 13},  // rank 7
        {12, 11}, // rank 8
        {12, 11}, // rank 9
        {-1, 14}, // rank 10
        {14, 18}, // rank 11
        {16, 15}, // rank 12
        {15, -1}, // rank 13
        {16, 17}, // rank 14
        {17, 18}, // rank 15
        {-1, -1}, // rank 16
        {-1, -1}, // rank 17
        {-1, -1}  // rank 18
    };

    int rank, num_proc; // process rank and total number of processes
    int low = -1;
    int high = -1; 
    int numbers[LEN_NUMBERS];

    MPI_Init(&argc, &argv); // initialize the MPI environment
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc); // get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get the rank of the process

    if (rank == ROOT) {
        //std::cout << "ROOT is going to distribute numbers now." << std::endl;

        std::ifstream numbers_file("numbers");
        if (!numbers_file.is_open() || numbers_file.bad()) {  // error opening numbers file
            std::cerr << "Could not open the file 'numbers'" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        //for (int i = 0; !numbers_file.eof(); i++) {
        for (int i = 0; i < LEN_NUMBERS; i++) {
            unsigned char number_c = numbers_file.get();
            int number_i = (int) number_c;
            int send_to_rank = i / 2;

            if (i != LEN_NUMBERS) {
                std::cout << number_i << " ";
            }
            else {
                std::cout << number_i;
            }
            MPI_Send(&number_i, 1, MPI_INT, send_to_rank, 0, MPI_COMM_WORLD); // Send number to correct process
        }
        std::cout << std::endl;
    }

    if (rank >= 0 && rank <= 3) { // TODO! better condition
        // Processes 0-4 will receive 2 numbers from rank 0
        //unsigned char tmp_low, tmp_high;
        //MPI_Request request1;

        MPI_Recv(&low, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //std::cout << "Rank " << rank << " received low=" << (int) low << " from " << ROOT << std::endl;
        MPI_Recv(&high, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //std::cout << "Rank " << rank << " received high=" << (int) high << " from " << ROOT << std::endl;

        // compare low and high
        if (low > high) {
            // exchange low and high
            int tmp = low;
            low = high;
            high = tmp;
        }

        // send low and high forward
        if (neighbour_matrix[rank][0] != -1) { // send low forward
            int send_low_to = neighbour_matrix[rank][0];
            MPI_Send(&low, 1, MPI_INT, send_low_to, 0, MPI_COMM_WORLD);
        }
        if (neighbour_matrix[rank][1] != -1) { // send low forward
            int send_high_to = neighbour_matrix[rank][0];
            MPI_Send(&high, 1, MPI_INT, send_high_to, 0, MPI_COMM_WORLD);
        }
    }

    //if (rank >= 4 && rank <= 9) { // TODO! better condition


    print_processor(rank, num_proc, low, high);

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}