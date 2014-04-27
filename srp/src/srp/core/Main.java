package srp.core;




/*
 * 
 -XX:CompileThreshold=50000 -XX:+CITime
  -XX:+DoEscapeAnalysis
-Djava.library.path=./lib -Xms256m -Xmx1024m

-XX:+PrintCompilation
-XX:+PrintCompilation -XX:+UnlockDiagnosticVMOptions -XX:+PrintInlining
*/



public class Main {

	/**
	 * @param args
	 * @throws Exception 
	 */
	static int RUNS = (int) 1e8;
	 private static final String FMT_DTL = "%-5s %6d %6.1f %6d %6.4f %6d %6.4f";
	 
	public static void main(String[] args) throws Exception {
		
		
		
		long time1 = System.currentTimeMillis();
		char[] cc = {'A','C'};
		char c = cc[0];
		for (int t = 0; t < 1e8; t++) {
//			BetaDistributionModel
		}
		long time2 = System.currentTimeMillis();

		System.out.println((time2 - time1) + "\t");
		
		 System.gc();
	        double x = 0;
	        long time = System.nanoTime();
	        for (int i = 0; i < RUNS; i++)
//	            x += StrictMath.log10(Math.PI + i/* 1.0 + i/1e9 */);
	        	x+=1;
	        long strictMath = System.nanoTime() - time;
	 
	        System.gc();
	        double y = 0;
	        time = System.nanoTime();
	        for (int i = 0; i < RUNS; i++)
	            y += Math.log(Math.PI + i/* 1.0 + i/1e9 */);
	        long fastTime = System.nanoTime() - time;
	 
	        System.gc();
	        double z = 0;
	        time = System.nanoTime();
	        for (int i = 0; i < RUNS; i++)
	            z += Math.log10(Math.PI + i/* 1.0 + i/1e9 */);
	        long mathTime = System.nanoTime() - time;
	 
	        report("log10", x + y + z, strictMath,fastTime,mathTime);	
	}
    private static void report(String name, double result, long strictMathTime, long fastMathTime, long mathTime) {
        long unitTime = strictMathTime;
        System.out.println("Keep HotSpot honest: " + result);
        System.out.println(String.format(FMT_DTL,
                name,
                strictMathTime / RUNS, (double) strictMathTime / unitTime,
                fastMathTime / RUNS, (double) fastMathTime / unitTime,
                mathTime / RUNS, (double) mathTime / unitTime
                ));
    }
}
