#include <stdio.h>
#include <time.h>
#include <stdlib.h>

/*with gcc and -O1 the recuresive func should look like:
00000000004006a2 <recursion> (File Offset: 0x6a2):
  4006a2:       test   %edi,%edi
  4006a4:       jne    4006a7 <recursion+0x5> (File Offset: 0x6a7)
  4006a6:       retq
  4006a7:       sub    $0x8,%rsp
  4006ab:       sub    $0x1,%edi
  4006ae:       callq  4006a2 <recursion> (File Offset: 0x6a2)
  4006b3:       add    $0x8,%rsp
  4006b7:       retq

hence, test+jne+sub(2x)+add+retq are llvm-estimated to be ~12 cycles
 -> rest should be callq and running it gives ~30 cycles per recursion
*/
void __attribute__ ((noinline)) recursion(int counter) {
    //printf("hello nr %d\n", counter);
    if (!counter)
        return;
    else
        recursion(counter-1); /* function calls itself */
}

int main(int argc, char **argv) {
    struct timespec tstart={0,0}, tend={0,0};
    clock_gettime(CLOCK_MONOTONIC, &tstart);
    recursion(atoi(argv[1]));
    clock_gettime(CLOCK_MONOTONIC, &tend);
    double t = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec)
        - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
    printf("time: %.8f s; at 2.5Ghz: avg. %.8f cycles per call\n",
            t, (2.5 * 1e9 * t) / atoi(argv[1]));
}
