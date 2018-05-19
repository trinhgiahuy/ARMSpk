#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <ittnotify.h>
#include <signal.h>
#include <stdlib.h>
#define STARTSDE(go,rank) {if(go && 0==rank && getenv("PCMPID")) kill(atoi(getenv("PCMPID")),SIGUSR1); __itt_resume(); __SSC_MARK(0x111);}
#define STOPSDE(go,rank) {__SSC_MARK(0x222); __itt_pause(); if(go && 0==rank && getenv("PCMPID")) kill(atoi(getenv("PCMPID")),SIGUSR1);}

#define RUM 10000000

int main() {
double mkrts, mkrte; // my kernel run-time
STOPSDE(0,0);

	float *a, *b;
	int i;
	MPI_Init(NULL, NULL);

	a=malloc(RUM * sizeof(*a));
	b=malloc(RUM  * sizeof(*b));
	for (i=0; i<RUM; a[i]=2, b[i]=3, i++);
	#pragma omp parallel for
	for (i=0; i<RUM; i++) a[i] *= b[i];
	for (i=0; i<RUM; b[i]=4, i++);

mkrts = MPI_Wtime();
STARTSDE(1,0);

	#pragma omp parallel for
	for (i=0; i<RUM; i++) a[i] *= b[i];

STOPSDE(1,0);
mkrte = MPI_Wtime();
printf("Walltime of the main kernel: %.6lf sec\n", mkrte - mkrts);

	printf("Hello, World %f %f\n", a[0], a[--i]);
	free(a);
	free(b);

	MPI_Finalize();

STARTSDE(0,0);
	return 0;
}
