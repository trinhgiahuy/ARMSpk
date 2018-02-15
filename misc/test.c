#include <stdio.h>
int main() {
	__SSC_MARK(0x222);
	float *a, *b;
	int i;
	a=malloc(1024 * sizeof(*a));
	b=malloc(1024 * sizeof(*b));
	for (i=0; i<1024; a[i]=2, b[i]=3, i++);
	for (i=0; i<1024; i++) a[i] *= b[i];
	for (i=0; i<1024; b[i]=4, i++);
	__SSC_MARK(0x111);
	for (i=0; i<1024; i++) a[i] *= b[i];
	__SSC_MARK(0x222);
	printf("Hello, World %f %f!", a[0], a[i]);
	free(a);
	free(b);
	__SSC_MARK(0x111);
	return 0;
}
