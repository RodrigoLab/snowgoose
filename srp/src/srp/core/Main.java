package srp.core;

import dr.inference.distribution.BetaDistributionModel;



/*
 -XX:CompileThreshold=50000 -XX:+CITime
  -XX:+DoEscapeAnalysis
-Djava.library.path=./lib -Xms256m -Xmx1024m
*/



public class Main {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		
		
		
		long time1 = System.currentTimeMillis();
		char[] cc = {'A','C'};
		char c = cc[0];
		for (int t = 0; t < 1e8; t++) {
//			BetaDistributionModel
		}
		long time2 = System.currentTimeMillis();

		System.out.println((time2 - time1) + "\t");

		
	}
}
