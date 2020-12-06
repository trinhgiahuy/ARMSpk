#include <mpi.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <stdlib.h>
//#define STARTSDE(go,rank) {__SSC_MARK(0x111);}
//#define STOPSDE(go,rank) {__SSC_MARK(0x222);}
//#include <ittnotify.h>
//#define STARTSDE(go,rank) {if(go && 0==rank && getenv("PCMPID")) kill(atoi(getenv("PCMPID")),SIGUSR1); __itt_resume(); __SSC_MARK(0x111);}
//#define STOPSDE(go,rank) {__SSC_MARK(0x222); __itt_pause(); if(go && 0==rank && getenv("PCMPID")) kill(atoi(getenv("PCMPID")),SIGUSR1);}

void __attribute__ ((noinline)) xxxx_add(int *i) {*i+=7;}

//https://crypto.stanford.edu/pbc/notes/pi/code.html
int __attribute__ ((noinline)) xxxx_compute(char *pi) {
    int r[2800 + 1];
    int i, j, k;
    int b, d;
    int c = 0, c0=(int)pi;

    for (i = 0; i < 2800; i++)
        r[i] = 2000;

    for (k = 2800, j = 0; k > 0; k -= 14, j+=4) {
        xxxx_add(&c0);
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
    return c0;
}

int __attribute__ ((noinline)) xxxx_sub(int c0) {
    int i, j, k;
    int d0 = -1 * c0;
    for (k = 2800, j = 0; k > 0; k -= 14, j+=4) {
        xxxx_add(&d0);
    }
    return -1 * d0;
}

int main(int argc, char** argv) {
    struct timespec tstart={0,0}, tend={0,0};
    clock_gettime(CLOCK_MONOTONIC, &tstart);
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //this 'first' region will not produce output, but lets just stop anyways
//    STOPSDE(0,world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    double t1 = MPI_Wtime();
//    STARTSDE(0,world_rank); //-> first region
    char* pi = malloc(801);
    pi[800] = '\0';
    int bs=xxxx_compute(pi);
    int bs2=xxxx_sub(bs);
    printf("compute Wtime=%e\n", MPI_Wtime()-t1);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors "
           "with 800 digits of PI: %s (bs:%d,%d)\n", processor_name, world_rank,
           world_size, pi, bs, bs2);
    printf("comp+print Wtime=%e\n", MPI_Wtime()-t1);

    t1 = MPI_Wtime();
    for (int i=0; i<1000000; i++)
        MPI_Bcast(pi, 800, MPI_CHAR, 0, MPI_COMM_WORLD);
//    STOPSDE(0,world_rank);
    printf("comm Wtime=%e (should be ~1.5s on 4 cores)\n", MPI_Wtime()-t1);

    free(pi);

//    STARTSDE(0,world_rank); //-> second region starts after here

    // Finalize the MPI environment.
    MPI_Finalize();
    clock_gettime(CLOCK_MONOTONIC, &tend);
    printf("whole app run took about %.8f seconds\n",
            ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) -
            ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
}
