#include <mpi.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <stdlib.h>
#define STARTSDE {__SSC_MARK(0x111);}
#define STOPSDE {__SSC_MARK(0x222);}

void __attribute__ ((noinline)) xxxx_add(int *i) {*i+=7;}
void __attribute__ ((noinline)) xxxx_sub(int *i) {*i-=7;}

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


    int pi=3-(7*3*2);
    STOPSDE; //ignore
    xxxx_add(&pi);
    xxxx_sub(&pi);
    STARTSDE;

    //count 5x add and 2x sub
    xxxx_add(&pi);
    xxxx_sub(&pi);
    xxxx_add(&pi);
    xxxx_add(&pi);
    xxxx_add(&pi);
    xxxx_add(&pi);
    xxxx_sub(&pi);

    STOPSDE;//ignore
    xxxx_add(&pi);
    xxxx_sub(&pi);
    STARTSDE;

    //count 5x add and 2x sub
    xxxx_add(&pi);
    xxxx_sub(&pi);
    xxxx_add(&pi);
    xxxx_add(&pi);
    xxxx_add(&pi);
    xxxx_add(&pi);
    xxxx_sub(&pi);

    STOPSDE;//ignore
    xxxx_add(&pi);
    xxxx_sub(&pi);
    STARTSDE;

    printf("rank %d: %d\n", world_rank, pi); //should print 3

    MPI_Finalize();
}
