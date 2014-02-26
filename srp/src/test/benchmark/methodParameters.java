package test.benchmark;

import java.util.Arrays;

import dr.math.MathUtils;

public class methodParameters {

	public static void main(String[] args) {
//		testPassing();
//		testComparing();
		testArraySwapVsCopy();
	}
	private static void testArraySwapVsCopy() {
		double zzz = 0;
		double[] a = new double[]{Math.random(), Math.random(), Math.random(), Math.random()}; 
		double[] b = new double[]{Math.random(), Math.random(), Math.random(), Math.random()};
		long time1 = System.currentTimeMillis();
		for (int t = 0; t < 1e9; t++) {
			double[] temp = a;
	    	a = b;
	    	b= temp;
	    	zzz+=1;
		}
		long time2 = System.currentTimeMillis();

		System.out.println((time2 - time1) + "\t");
		System.out.println(Arrays.toString(a) +"\t"+ Arrays.toString(b));
		
		time1 = System.currentTimeMillis();
		for (int t = 0; t < 1e9; t++) {
	    	System.arraycopy(a, 0, b, 0, 4);
	    	zzz+=1;
		}
		time2 = System.currentTimeMillis();

		System.out.println((time2 - time1) + "\t");
		System.out.println(zzz);
		time1 = System.currentTimeMillis();
		for (int t = 0; t < 1e9; t++) {
			zzz += 1;
		}
		time2 = System.currentTimeMillis();

		System.out.println((time2 - time1) + "\t");
		System.out.println(zzz);
	}
	private static void testComparing() {
		int ite =  (int) 1e9;
		int result = 0;
		double tt = 0;
		long time1 = System.currentTimeMillis();
		for (int t = 0; t < ite; t++) {
			tt = t;
			if(t != 100){
				result += tt;
			}
		}
		System.out.println(tt);
		long time2 = System.currentTimeMillis();

		System.out.println("Time: "+(time2 - time1) + "\t" + (time2 - time1)/ite +"/caluclation");
		System.out.println(result);
		result = 0;
		time1 = System.currentTimeMillis();
		for (int t = 0; t > -ite; t--) {
			tt = t;
			if(tt != 100000.1){
				result += tt;
			}
		}
		time2 = System.currentTimeMillis();
		System.out.println(result);
		System.out.println("Time: "+(time2 - time1) + "\t" + (time2 - time1)/ite +"/caluclation");
		
	}
	private static void testPassing(){
		double result = 0;
		int[] indexArray = new int[200];
		for (int i = 0; i < indexArray.length; i++) {
			indexArray[i] = i;//MathUtils.nextInt(100);
		}
		double[] temp = new double[100];
		for (int i = 0; i < temp.length; i++) {
			temp[i] = Math.random();
		}
		int ite = (int) 1e8;
		int index = 0;
		long time1 = System.currentTimeMillis();
		
			//		for (int tt = 0; tt < 10; tt++) {
			for (int t = 0; t < ite; t++) {
//				result = passArray(indexArray, index, temp);
				result = passArray(indexArray);
			}
//		}
		long time2 = System.currentTimeMillis();
		System.out.println(result);
		System.out.println("Time: "+(time2 - time1) + "\t");
		
		time1 = System.currentTimeMillis();
//		for (int tt = 0; tt < 10; tt++) {
			for (int t = 0; t < ite; t++) {
//				int x = indexArray[index];
//				result = passIndex(x, temp);
//				result = passIndex(temp[x]);
				result = passInt(0, indexArray.length);
			}
//		}
		time2 = System.currentTimeMillis();
		System.out.println(result);
		System.out.println("Time: "+(time2 - time1) + "\t");
		
		
		
	}
	
	private static double passIndex(double d) {
		return d+1;
	}

	private static double passArray(int[] indexArray, int index, double[] temp) {
		int x = indexArray[index];
		return temp[x]+1;
	}
	
	private static double passIndex(int index, double[] temp) {
		return temp[index]+1;
	}
	
	private static int passInt(int start, int length){
		int k= 0;
		for (int s = 0; s < length; s++) {
			k += start+s;
		}
		return k;
	}
	
	private static int passArray(int[] indexArray){
		int k= 0;
		for (int i = 0; i < indexArray.length; i++) {
			k += indexArray[i];
		}
		return k;
	}
}
