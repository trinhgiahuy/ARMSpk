#include <stdio.h>
#include <mpi.h>

#define STARTSDE __SSC_MARK(0x111);
#define STOPSDE __SSC_MARK(0x222);

int main() {
double mkrts, mkrte; // my kernel run-time
STOPSDE;

	float *a, *b;
	int i;
	MPI_Init(NULL, NULL);

	a=malloc(10240 * sizeof(*a));
	b=malloc(10240 * sizeof(*b));
	for (i=0; i<10240; a[i]=2, b[i]=3, i++);
	#pragma omp parallel for
	for (i=0; i<10240; i++) a[i] *= b[i];
	for (i=0; i<10240; b[i]=4, i++);

mkrts = MPI_Wtime();
STARTSDE;

	#pragma omp parallel for
	for (i=0; i<10240; i++) a[i] *= b[i];

STOPSDE;
mkrte = MPI_Wtime();
printf("Walltime of the main kernel: %.6lf sec\n", mkrte - mkrts);

	printf("Hello, World %f %f\n", a[0], a[--i]);
	free(a);
	free(b);

	MPI_Finalize();

STARTSDE;
	return 0;
}
