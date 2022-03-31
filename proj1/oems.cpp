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
        std::ifstream numbers_file("numbers");
        if (!numbers_file.is_open() || numbers_file.bad()) {  // error opening numbers file
            std::cerr << "Could not open the file 'numbers'" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

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

    MPI_Recv(&low, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&high, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

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
        //std::cout << "Sending from " << rank << " low to " << send_low_to << std::endl;
        MPI_Send(&low, 1, MPI_INT, send_low_to, 0, MPI_COMM_WORLD);
    }

    if (neighbour_matrix[rank][1] != -1) { // send low forward
        int send_high_to = neighbour_matrix[rank][1];
        //std::cout << "Sending from " << rank << " high to " << send_high_to << std::endl;
        MPI_Send(&high, 1, MPI_INT, send_high_to, 0, MPI_COMM_WORLD);
    }

    //print_processor(rank, num_proc, low, high);

    switch (rank) {
        case 10: 
            MPI_Send(&low, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            break;
        case 16:
        case 17:
        case 18: 
            MPI_Send(&low, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&high, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            break;
        case 13:
            MPI_Send(&high, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            break;
    }

    if (rank == 0) {
        for (int i = 0; i < num_proc; i++) {
            switch (i) {
                case 10: 
                    MPI_Recv(&numbers[0], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                    break;
                case 16:
                    MPI_Recv(&numbers[1], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                    MPI_Recv(&numbers[2], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    break; 
                case 17:
                    MPI_Recv(&numbers[3], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                    MPI_Recv(&numbers[4], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    break; 
                case 18: 
                    MPI_Recv(&numbers[5], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                    MPI_Recv(&numbers[6], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    break; 
                case 13:
                    MPI_Recv(&numbers[7], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                    break;
            }
        }

        for (int i = 0; i < LEN_NUMBERS; i++) {
            std::cout << numbers[i] << std::endl;
        }
    }


    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}