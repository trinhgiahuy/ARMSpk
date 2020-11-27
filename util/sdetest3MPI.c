#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    for (int i=0; i<10; i++)
        MPI_Bcast(&world_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Print off a hello world message
    printf("Hello world from rank %d out of %d\n", world_rank, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();
}
