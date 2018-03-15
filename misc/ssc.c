#include <ittnotify.h>
#define STARTSDE {__itt_resume(); __SSC_MARK(0x111);}
#define STOPSDE {__itt_pause(); __SSC_MARK(0x222);}

int ssc_mark_start (int x) {
	STARTSDE;
	return (x << 0x1);
}
int ssc_mark_stop (int y) {
	STOPSDE;
	return (y >> 0x2);
}
