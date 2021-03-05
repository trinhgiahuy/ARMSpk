#include <SdeTestJava.h>

#include <ittnotify.h>
#include <signal.h>
#include <stdlib.h>
#define STARTSDE(go,rank) {if(go && 0==rank && getenv("PCMPID")) kill(atoi(getenv("PCMPID")),SIGUSR1); __itt_resume(); __SSC_MARK(0x111);}
#define STOPSDE(go,rank) {__SSC_MARK(0x222); __itt_pause(); if(go && 0==rank && getenv("PCMPID")) kill(atoi(getenv("PCMPID")),SIGUSR1);}

JNIEXPORT jint JNICALL
Java_SdeTestJava_ssc_1mark_1start(JNIEnv *env, jobject obj, jint go, jint rank)
{
	STARTSDE(go,rank);
	return (rank << 0x1);
}
JNIEXPORT jint JNICALL
Java_SdeTestJava_ssc_1mark_1stop(JNIEnv *env, jobject obj, jint go, jint rank)
{
	STOPSDE(go,rank);
	return (rank >> 0x2);
}
