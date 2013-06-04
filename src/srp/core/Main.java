package srp.core;
/*
 -XX:CompileThreshold=50000 -XX:+CITime
  -XX:+DoEscapeAnalysis
-Djava.library.path=./lib -Xms256m -Xmx1024m
*/

import java.util.Arrays;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.HaplotypeModelUtils;
import test.mcmc.MCMCTrueTree;
import dr.evolution.alignment.Alignment;
import dr.math.MathUtils;


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
//			int a = MathUtils.nextInt(); //8-9k
//			if(a>0){
//				c = cc[0];
//			}
//			else{
//				c = cc[1];
//			}
//			boolean b = MathUtils.nextBoolean(); //9k
//			if(b){
//				c =cc[0];
//			}
//			else{
//				c = cc[1];
//			}

			//			MathUtils.nextDouble(); //10k
//			MathUtils.nextByte(); //9k
			
//			int a = MathUtils.nextInt(2); //9.5k
//			c = cc[a];
			
//			System.out.println(MathUtils.nextByte());
		}
		long time2 = System.currentTimeMillis();

		System.out.println((time2 - time1) + "\t");

		
	}
}
