int ssc_mark_start (int x) {
	__SSC_MARK(0x111);
	return (x << 0x1);
}
int ssc_mark_stop (int y) {
	__SSC_MARK(0x222);
	return (y >> 0x2);
}
