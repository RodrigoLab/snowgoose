package srp.likelihood.haplotypes;


public class LikelihoodUtils {

	private static final int ZERO = 0;
	private static final int ONE = 1;


	public static int hamDist(String s1, String s2)  {
//
//		if (s1.length() != s2.length()) {
//			System.err.println("Different length\t"+  s1.length() +"\t"+ s2.length()) ;
//			System.exit(-1);
//		}
		 
		int count = 0;
		int a;
		for (int i = 0; i < s1.length(); i++) {

			a = (s1.charAt(i) - s2.charAt(i));
			count+= (a==ZERO)? 0:1;
		}
		return count;
	}



	public static int hamDist(char[] c1, char[] c2) {


		int count = 0;

		for (int i = 0; i < c1.length; i++) {
//			System.out.println(c1[i] +"\t"+ c2[i] +"\t"+ (c1[i]-c2[i]));
			if ( (c1[i]-c2[i])!=0 ){
//			if (s1.charAt(i) != s2.charAt(i) ){
				count++;
			}
		}
		return count;
	}





	static final Byte bb = new Byte((byte) 0);
	
	
	public static int hamDist(char[] srCharArray, char[] hapCharArray, int start) {
		
		int count = 0;
		int a;
		
		for (int i = 0; i < srCharArray.length; i++) {
			a = (srCharArray[i] - hapCharArray[i+start]);
			
			count += (a==ZERO)? 0:1;

		}

		return count;
	}

//static int cc = 500;
//static int c2 = 0;//
	static int counts[] = new int[1000];
	static char[] srCharArray;
	static char[] hapCharArray;
	static int a, noComb;
	static char cc;
	private static String hapString;
	private static String srpString;
	
	static int match;
	
	public static int[] hamDistAll(char[] srCharArray2, char[] hapCharArray2) {
		srCharArray = srCharArray2;
		hapCharArray = hapCharArray2;
		
		int noComb = hapCharArray.length - srCharArray.length + 1; 		
//				int noComb = hLength - srLength + 1;
//System.out.println(noComb);
//		Arrays.fill(counts, 0);
//		int x = noComb*100;
//		int counts[] = new int[550];
//		if (noComb>counts.length) {
//			Arrays.fill(counts, 1);
//			c2++;
//			return counts;
//		counts = new int[noComb];
//			Arrays.fill(counts, 0);
//		int[] counts = Arrays.copyOfRange(counts2, 0, noComb);
//		int[] counts = Arrays.copyOf(counts2, noComb);
//		Arrays.fill(counts, 0, noComb, 0);
//		int[] counts = ArrayUtils.clone(counts2);
//		Ints.
		
//			int[] x = ArrayUtils.EMPTY_INT_ARRAY;
//			Ints.ensureCapacity(array, minLength, padding)
//			System.out.println(x.length);
//			System.out.println(noComb);
//		}
//		else {
//			c2++;
//		int a;
//		char cc;
			for (int i = 0; i < noComb; i++) {
				counts[i]=0;
			}
//		}

		
		for (int j = 0; j < srCharArray.length; j++) {
			cc = srCharArray[j];
			for (int i = 0; i < noComb; i++) {
				int a = (cc - hapCharArray[j+i]);
				counts[i] += (a==ZERO)? ZERO:ONE;
			}
		}

		return counts;
	}



	public static int[] hamDistAll(char[] srCharArray2, String hap) {
		srCharArray = srCharArray2;
		hapString = hap;
		
		int noComb = hapString.length() - srCharArray.length + 1; 		
//				int noComb = hLength - srLength + 1;
//System.out.println(noComb);
//		Arrays.fill(counts, 0);
//		int x = noComb*100;
//		int counts[] = new int[550];
//		if (noComb>counts.length) {
//			Arrays.fill(counts, 1);
//			c2++;
//			return counts;
//		counts = new int[noComb];
//			Arrays.fill(counts, 0);
//		int[] counts = Arrays.copyOfRange(counts2, 0, noComb);
//		int[] counts = Arrays.copyOf(counts2, noComb);
//		Arrays.fill(counts, 0, noComb, 0);
//		int[] counts = ArrayUtils.clone(counts2);
//		Ints.
		
//			int[] x = ArrayUtils.EMPTY_INT_ARRAY;
//			Ints.ensureCapacity(array, minLength, padding)
//			System.out.println(x.length);
//			System.out.println(noComb);
//		}
//		else {
//			c2++;
//		int a;
//		char cc;
			for (int i = 0; i < noComb; i++) {
				counts[i]=0;
			}
//		}

		
		for (int j = 0; j < srCharArray.length; j++) {
			cc = srCharArray[j];
			for (int i = 0; i < noComb; i++) {
				int a = (cc - hapString.charAt(j+i));
				counts[i] += (a==ZERO)? ZERO:ONE;
			}
		}

		return counts;

	}

	public static int Dist(int start, int end, String srp, String hap) {
//		hapCharArray = hapString.toCharArray();
		hapString = hap;
		srpString = srp;
		
		int count = 0;
		int dist = 0;
//		System.out.println(Arrays.toString(srCharArray));
//		char[] c = 		ArrayUtils.subarray(hapCharArray, start, end);
//		System.out.println(Arrays.toString(c));
		for (int i = start; i < end; i++) {
//			match = (srCharArray[count] - hapCharArray[i]);
			match = (srpString.charAt(count) - hapString.charAt(i));
			dist += (match==ZERO)? 0:1;
//				System.out.println(count +"\t"+ i+"\t"+ match +"\t"+ (srCharArray[count] - hapCharArray[i]) +"\t");
			count++;
		}
//		System.out.println(dist);
//		System.out.println("===============");
		return dist;
	}
	public static int Dist(int start, int end, char[] srp, String hap) {
//		hapCharArray = hapString.toCharArray();
		hapString = hap;
		srCharArray = srp;
		
		int count = 0;
		int dist = 0;
//		System.out.println(Arrays.toString(srCharArray));
//		char[] c = 		ArrayUtils.subarray(hapCharArray, start, end);
//		System.out.println(Arrays.toString(c));
		for (int i = start; i < end; i++) {
//			match = (srCharArray[count] - hapCharArray[i]);
			match = (srCharArray[count] - hapString.charAt(i));
			dist += (match==ZERO)? 0:1;
//				System.out.println(count +"\t"+ i+"\t"+ match +"\t"+ (srCharArray[count] - hapCharArray[i]) +"\t");
			count++;
		}
//		System.out.println(dist);
//		System.out.println("===============");
		return dist;
	}
	public static int Dist(int start, int end, char[] srp, char[] hapChar) {
		hapCharArray = hapChar;
		srCharArray = srp;
		
		int count = 0;
		int dist = 0;
//		System.out.println(Arrays.toString(srCharArray));
//		char[] c = 		ArrayUtils.subarray(hapCharArray, start, end);
//		System.out.println(Arrays.toString(c));
		for (int i = start; i < end; i++) {
			match = (srCharArray[count] - hapCharArray[i]);
			dist += (match==ZERO)? 0:1;
//				System.out.println(count +"\t"+ i+"\t"+ match +"\t"+ (srCharArray[count] - hapCharArray[i]) +"\t");
			count++;
		}
//		System.out.println(dist);
//		System.out.println("===============");
		return dist;
	}



}
