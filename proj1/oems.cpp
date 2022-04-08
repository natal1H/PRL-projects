/**
 * @file oems.cpp
 * @author Natália Holková (xholko02@stud.fit.vutbr.cz)
 * @brief PRL - project 1 - Odd-even merge sort for 8 numbers
 * @date 2022-04-08 
 */

#include <iostream>
#include <fstream>
#include <mpi.h>

#define ROOT 0
#define LEN_NUMBERS 8
#define INPUT_FILENAME "numbers"

using namespace std;

/**
 * @brief Reads 8 1 byte numbers from INPUT_FILENAME, prints them and sends directly to processors. 
 */
void distribute_numbers() {
    ifstream numbers_file(INPUT_FILENAME);
    if (!numbers_file.is_open() || numbers_file.bad()) {  // Error opening numbers file
        cerr << "Could not open the file 'numbers'" << endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int i = 0; i < LEN_NUMBERS; i++) {
        unsigned char number_c = numbers_file.get();
        int number_i = (int) number_c;
        int send_to_rank = i / 2;

        if (i != LEN_NUMBERS)
            cout << number_i << " ";
        else 
            cout << number_i;


        MPI_Send(&number_i, 1, MPI_INT, send_to_rank, 0, MPI_COMM_WORLD); // Send number to correct process
    }
    cout << endl;
}

/**
 * @brief Procesor receives two values and swaps them if necessary
 * 
 * @param low_p Pointer to where low value is stored
 * @param high_p Pointer to where high value is stored
 */
void receive_inputs(int *low_p, int *high_p) {
    MPI_Recv(low_p, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(high_p, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // compare low and high and exchange 
    if (*low_p > *high_p)
        swap(*low_p, *high_p);
}

/**
 * @brief Processor send provided value to specified rank
 * 
 * @param to_rank Rank of receiver
 * @param val Value to send forward
 */
void send_value_forward(int to_rank, int val) {
    MPI_Send(&val, 1, MPI_INT, to_rank, 0, MPI_COMM_WORLD);
}

/**
 * @brief Processor goes through all processors and waits to receive message from some of them. 
 *        Only ROOT processor should call this function.
 * 
 * @param arr Pointer to array where values will be stored
 * @param num_proc Total number of processors
 */
void collect_sorted_values(int *arr, int num_proc) {
    for (int i = 0; i < num_proc; i++) {
        switch (i) {
            case 10: 
                MPI_Recv(&arr[0], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                break;
            case 16:
                MPI_Recv(&arr[1], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Recv(&arr[2], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                break; 
            case 17:
                MPI_Recv(&arr[3], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Recv(&arr[4], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                break; 
            case 18: 
                MPI_Recv(&arr[5], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Recv(&arr[6], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                break; 
            case 13:
                MPI_Recv(&arr[7], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                break;
        }
    }
}

// MAIN PROGRAM

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
        {0, 14}, // rank 10
        {14, 18}, // rank 11
        {16, 15}, // rank 12
        {15, 0}, // rank 13
        {16, 17}, // rank 14
        {17, 18}, // rank 15
        {0, 0}, // rank 16
        {0, 0}, // rank 17
        {0, 0}  // rank 18
    };

    int rank, num_proc; // process rank and total number of processes
    int low = -1;
    int high = -1; 

    MPI_Init(&argc, &argv); // initialize the MPI environment
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc); // get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get the rank of the process

    // ROOT reads and distributes input
    if (rank == ROOT)
        distribute_numbers();

    // Each processor receives two inputs
    receive_inputs(&low, &high);

    // send low and high forward
    send_value_forward(neighbour_matrix[rank][0], low);
    send_value_forward(neighbour_matrix[rank][1], high);

    // ROOT will collect sorted numbers and display
    if (rank == ROOT) {
        int numbers[LEN_NUMBERS] = {-1};
        collect_sorted_values(numbers, num_proc);
        
        for (int i = 0; i < LEN_NUMBERS; i++)
            cout << numbers[i] << endl;
    }

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}