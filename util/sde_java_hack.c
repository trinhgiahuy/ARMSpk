//icc -Wall -fPIC -c handler.c ; ld -shared handler.o -o handler.so
//LD_PRELOAD=handler.so java
//XXX: don't use global variables in here!!!

#include <signal.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

#include <ittnotify.h>
#include <stdlib.h>
#define STARTSDE(go,rank) {if(go && 0==rank && getenv("PCMPID")) kill(atoi(getenv("PCMPID")),SIGUSR1); __itt_resume(); __SSC_MARK(0x111);}
#define STOPSDE(go,rank) {__SSC_MARK(0x222); __itt_pause(); if(go && 0==rank && getenv("PCMPID")) kill(atoi(getenv("PCMPID")),SIGUSR1);}

void start_sde(int sig) {
    FILE *log = fopen("/dev/shm/sde.java.hack.log", "a");
    assert(log);
    fprintf(log, "starting @%f\n", (double)clock() / CLOCKS_PER_SEC);
    fclose(log);
    STARTSDE(0,0);
}

void stop_sde(int sig) {
    STOPSDE(0,0);
    FILE *log = fopen("/dev/shm/sde.java.hack.log", "a");
    assert(log);
    fprintf(log, "stopping @%f\n", (double)clock() / CLOCKS_PER_SEC);
    fclose(log);
}

void _init(void)
{
    remove("/dev/shm/sde.java.hack.log");
    FILE *pid = fopen("/dev/shm/sde.java.hack.pid", "w");
    assert(pid);
    fprintf(pid, "%lu\n", (unsigned long)getpid());
    fclose(pid);
    signal(SIGUSR1, start_sde);
    signal(SIGUSR2, stop_sde);
}
