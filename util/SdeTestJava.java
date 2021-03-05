public class SdeTestJava {

	//============ C interface =========================
	public native int ssc_mark_start (int go, int rank);
	public native int ssc_mark_stop (int go, int rank);
	static {
		System.loadLibrary("ssc");
	}
	//==================================================

	private void say_hello() {
		ssc_mark_start(0,0);
		System.out.println("Hello World");
		ssc_mark_stop(0,0);
	}

	public static void main(String[] args) {
		new SdeTestJava().say_hello();
	}
}

