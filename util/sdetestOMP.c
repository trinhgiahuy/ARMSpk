#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef ARRAY_SIZE
#define ARRAY_SIZE 1024*1024*1024
#endif

int main(int argc, char** argv) {
    long i = 0;
    double *a = malloc(ARRAY_SIZE * sizeof(*a));

    #pragma omp barrier
    #pragma omp parallel for schedule(static)
    for (i = 0; i < ARRAY_SIZE; i++)
        a[i] = sqrt(cos(i)) / sqrt(sin(i));
    #pragma omp barrier

    #pragma omp parallel
    printf("hello world from thread %d\n", omp_get_thread_num());
    free(a);

    return 0;
}
