#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#ifdef __GNUC__
static inline void *__movsb(void *d, const void *s, size_t n)
{
    asm("rep movsb"
            : "=D" (d),
            "=S" (s),
            "=c" (n)
            : "0" (d),
            "1" (s),
            "2" (n)
            : "memory");
    return d;
}
#else
# error "no compat helper for your compiler"
#endif

int main(int argc, char **argv) {
    struct timespec tstart={0,0}, tend={0,0};
    int len = atoi(argv[1]);
    char *s1 = calloc(len, 1);
    char *s2 = calloc(len, 1);
    memcpy(s2, argv[2], len);
    clock_gettime(CLOCK_MONOTONIC, &tstart);
    for (int i=0; i<1e6; i++) __movsb(s1, s2, len);
    clock_gettime(CLOCK_MONOTONIC, &tend);
    double t = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec)
               - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);

    printf("%s (time: %.8f s; at 2.5Ghz: avg. %.8f cycles per movsb)\n",
           s1, t, (2.5 * 1e9 * t) / (len * 1e6));
}
