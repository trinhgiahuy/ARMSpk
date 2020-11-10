#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

//https://crypto.stanford.edu/pbc/notes/pi/code.html
void compute(char *pi) {
    int r[2800 + 1];
    int i, j, k;
    int b, d;
    int c = 0;

    for (i = 0; i < 2800; i++)
        r[i] = 2000;

    for (k = 2800, j = 0; k > 0; k -= 14, j+=4) {
        d = 0;

        i = k;
        for (;;) {
            d += r[i] * 10000;
            b = 2 * i - 1;

            r[i] = d % b;
            d /= b;
            i--;
            if (i == 0) break;
            d *= i;
        }
        sprintf(pi+j, "%04d", c + d / 10000);
        c = d % 10000;
    }
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

    char* pi = malloc(801);
    pi[800] = '\0';
    compute(pi);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors "
           "with 800 digits of PI: %s\n", processor_name, world_rank,
           world_size, pi);

    free(pi);

    // Finalize the MPI environment.
    MPI_Finalize();
}
