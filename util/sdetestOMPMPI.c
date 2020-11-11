#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define STREAM_ARRAY_SIZE 1024*1024

void __attribute__ ((noinline)) xxxx_init(double *a, double *b, double *c) {
    int j;
    #pragma omp parallel for
    for (j=0; j<STREAM_ARRAY_SIZE; j++) {
        a[j] = 1.0;
        b[j] = 2.0;
        c[j] = 0.0;
    }
}

void __attribute__ ((noinline)) xxxx_stream(double *a, double *b, double *c,
                                            double scalar) {
    int j;
    #pragma omp parallel for
    for (j = 0; j < STREAM_ARRAY_SIZE; j++)
        c[j] = a[j] + scalar * b[j];
}

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    double *a = malloc(STREAM_ARRAY_SIZE * sizeof(*a));
    double *b = malloc(STREAM_ARRAY_SIZE * sizeof(*b));
    double *c = malloc(STREAM_ARRAY_SIZE * sizeof(*c));

    xxxx_init(a, b, c);
    xxxx_stream(a, b, c, (double)(world_rank+1));

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors "
           "with first/last values of stream array C: %f %f\n", processor_name,
           world_rank, world_size, c[0], c[STREAM_ARRAY_SIZE-1]);

    free(a); free(b); free(c);

    // Finalize the MPI environment.
    MPI_Finalize();
}
