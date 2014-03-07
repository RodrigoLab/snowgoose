package test.benchmark;

import java.util.Arrays;
import java.util.BitSet;
import java.util.Random;

import cern.colt.bitvector.BitVector;
import dr.math.MathUtils;

public class methodParameters {

	public static void main(String[] args) {
//		testPassing();
//		testComparing();
//		testArraySwapVsCopy();
		
//		int[] test = new int[]{0,1,2,3,4,5};
//		ShuffleArray(test);
//		System.out.println(Arrays.toString(test));
		
		testBit();
	}
	
	private static void testBit(){
		
		BitSet bs = new BitSet();
		bs.flip(10);
		bs.flip(2);
		System.out.println(bs.cardinality());
		System.out.println(bs.toString());
		System.out.println(Arrays.toString(bs.toLongArray()));
		System.out.println(bs.cardinality());
		bs.or( bs.valueOf(new long[]{512}) );
		System.out.println(bs.toString());
		System.out.println(Arrays.toString(bs.toLongArray()));
		bs.set(2, 5, true);
		System.out.println(bs.toString());
		bs.set(3, 4, false);
		System.out.println(bs.toString());
		System.out.println(bs.nextClearBit(2));
		System.out.println(bs.nextSetBit(2));
		
		bs.set(20);
		System.out.println(bs.toString());
		System.out.println(Arrays.toString(bs.toLongArray()));
		System.out.println(bs.nextSetBit(10));
		System.out.println(bs.nextSetBit(11));
		bs.set(70);
		System.out.println(bs.toString());
		System.out.println(Arrays.toString(bs.toLongArray()));
		bs.clear();
		bs.set(0);
		System.out.println(bs.toString());
		System.out.println(Arrays.toString(bs.toLongArray()));
		bs.set(1);
		System.out.println(bs.toString());
		System.out.println(Arrays.toString(bs.toLongArray()));
		bs.set(100);
		int index = -1;
		do{
			index = bs.nextSetBit(index+1);
			System.out.println(index);
		}while(index!= -1);
		
	
		System.out.println("===========");
		
		BitVector bitV = new BitVector(11);
		bitV.set(0);
		System.out.println(bitV.toString());
		System.out.println(bitV.getLongFromTo(0, 10));
		bitV.set(1);
		System.out.println(bitV.toString());
		System.out.println(bitV.getLongFromTo(0, 10));
		bitV.clear();
		bitV.set(0);
		bitV.set(3);
		bitV.set(5);
		bitV.set(9);
		System.out.println(bitV.toString());
		System.out.println(bitV.getLongFromTo(0, 10));
		System.out.println(bitV.size());
		int index2 = -1;
		do{
			index2 +=1;
			index2 = bitV.indexOfFromTo(index2, 10, true);
			System.out.println(index2);
//			index +=1;
		}while(index2 >=0 && index2<10);
	}
	
	private static void ShuffleArray(int[] array)
	{
	    int index;
	    Random random = new Random();
	    for (int i = array.length - 1; i > 0; i--)
	    {
	        index = random.nextInt(i + 1);
	        if (index != i)
	        {
	        	System.out.println(array[index] +"\t"+  array[i]);
	            array[index] ^= array[i];
	            System.out.println(array[index] +"\t"+  array[i]);
	            array[i] ^= array[index];
	            System.out.println(array[index] +"\t"+  array[i]);
	            array[index] ^= array[i];
	            System.out.println(array[index] +"\t"+  array[i]);
	            System.out.println();
	        }
	    }
//	    for (int i = array.length - 1; i > 0; i--)
//	    {
//	        index = random.nextInt(i + 1);
//	        temp = array[index];
//	        array[index] = array[i];
//	        array[i] = temp;
//	    }
	}
	private static void FY_insideOut(int n, int k){
		
////		FY shuffle
//		To shuffle an array a of n elements (indices 0..n-1):
		  
//			The "inside-out" algorithm[edit]
//
//			The Fisher–Yates shuffle, as implemented by Durstenfeld, is an in-place shuffle. That is, given a preinitialized array, it shuffles the elements of the array in place, rather than producing a shuffled copy of the array. This can be an advantage if the array to be shuffled is large.
//
//			To simultaneously initialize and shuffle an array, a bit more efficiency can be attained by doing an "inside-out" version of the shuffle. In this version, one successively places element number i into a random position among the first i positions in the array, after moving the element previously occupying that position to position i. In case the random position happens to be number i, this "move" (to the same place) involves an uninitialised value, but that does not matter, as the value is then immediately overwritten. No separate initialization is needed, and no exchange is performed. In the common case where source is defined by some simple function, such as the integers from 0 to n - 1, source can simply be replaced with the function since source is never altered during execution.
//
//			To initialize an array a of n elements to a randomly shuffled copy of source, both 0-based:
//		math.su
		int[] source = new int[]{0,1,2,3,4,5};
		k = 5;	
		int[] a = new int[k];
		for (int i = n-1; i >0; i--) {
			int j = MathUtils.nextInt(n);
			a[k] = source[j];
			source[j] = source[i];
		}
//		for i from n − 1 downto 1 do
//		       j ← random integer with 0 ≤ j ≤ i
//		       exchange a[j] and a[i]
		    	
		
		
		a[0] = source[0];
		for (int i = 0; i < a.length; i++) {
			int j = MathUtils.nextInt(n);
			if (j != i)
	          a[i] = a[j];
			a[j] = source[i];
			n--;
		}
//			  for i from 1 to n − 1 do
//			      j ← random integer with 0 ≤ j ≤ i
//			      if j ≠ i
//			          a[i] ← a[j]
//			      a[j] ← source[i]
//			    		  http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
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
