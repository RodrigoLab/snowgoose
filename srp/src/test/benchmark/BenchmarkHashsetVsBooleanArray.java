package test.benchmark;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import srp.core.DataImporter;
import srp.shortreads.AlignmentMapping;

import com.google.common.primitives.Ints;

import dr.evolution.alignment.Alignment;
import dr.math.MathUtils;

/*
 * 
 * count: number of random numbers, depth for each position
 * Element: total size of array
For count: 100		Element:100		Use HashSet:	152		0.0152
For count: 1000		Element:1000	Use HashSet:	698		0.0698
For count: 1000		Element:10000	Use HashSet:	1255	0.1255
For count: 1000		Element:100000	Use HashSet:	947		0.0947
For count: 10000	Element:1000000	Use HashSet:	9213	0.9213
For count: 100000	Element:10000	Use HashSet:	54520	5.452



For count: 100		Element:100			Use boolean[]:	99	0.0099
For count: 1000		Element:1000		Use boolean[]:	275	0.0275
For count: 1000		Element:10000		Use boolean[]:	444	0.0444
For count: 1000		Element:100000		Use boolean[]:	624	0.0624

For count: 10000	Element:10000		Use boolean[]:	2965	0.2965
For count: 10000	Element:100000		Use boolean[]:	3163	0.3163
For count: 10000	Element:1000000		Use boolean[]:	6523	0.6523
For count: 10000	Element:10000000	Use boolean[]:	49110	4.911

For count: 100000	Element:1000		Use boolean[]:	27673	2.7673
For count: 100000	Element:10000		Use boolean[]:	27644	2.7644
For count: 100000	Element:100000		Use boolean[]:	28021	2.8021


For count: 100		Element:100			Reuse HashSet:	170	0.017
For count: 1000		Element:1000		Reuse HashSet:	697	0.0697
For count: 1000		Element:10000		Reuse HashSet:	679	0.0679
For count: 1000		Element:100000		Reuse HashSet:	669	0.0669
For count: 1000		Element:1000000		Reuse HashSet:	892	0.0892
For count: 1000		Element:10000000	Reuse HashSet:	858	0.0858
For count: 1000		Element:10000000	Reuse HashSet:	858	0.0858

For count: 10000	Element:10000		Reuse HashSet:	6210	0.621
For count: 10000	Element:100000		Reuse HashSet:	7051	0.7051
For count: 10000	Element:1000000		Reuse HashSet:	6757	0.6757
For count: 10000	Element:10000000	Reuse HashSet:	6550	0.655

For count: 100000	Element:1000		Reuse HashSet:	37950	3.795
For count: 100000	Element:10000		Reuse HashSet:	47168	4.7168
For count: 100000	Element:100000		Reuse HashSet:	97540	9.754

700 reads, ~130-150x in the middle
AlignmentMapping with boolean[]:	3274	0.003274
AlignmentMapping with HashMap:		30020	0.03002

5000 reads
AlignmentMapping with boolean[]:	3540	0.0354
AlignmentMapping with HashMap:	26402	0.26402

1M reads??
*/


public class BenchmarkHashsetVsBooleanArray {

	public static void main(String[] args) throws Exception {

//		int ite = (int) 1e5;
		basicTest();
//		
//		AlignmentMapping aMap = setup();
//		useBooleanMapping(ite, 12, aMap);
//		useHashMapMapping(ite, 12, aMap);
	}
	
	private static void useBooleanMapping(int ite, int swapLength, AlignmentMapping aMap) {

//		System.out.print("For count: "+totalCount +"\tElement:"+ totalElements +"\t");
		
		long time1 = System.currentTimeMillis();
		boolean[] srpSwitch = new boolean[aMap.getSrpCount()];
		int[] twoPositions = new int[2];
		int hapLength = aMap.getLength();
		for (int t = 0; t < ite; t++) {
			twoPositions[0] = MathUtils.nextInt(hapLength - swapLength);
			twoPositions[1]= twoPositions[0] + swapLength;
			
			Arrays.fill(srpSwitch, false);
			for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
				
				ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);
//				System.out.println("Site: "+k +"\t"+ mapToSrp.size());
				for (int i : mapToSrp) {
					srpSwitch[i] = true;
				}
			}
			
//			for (int i = 0; i < totalCount; i++) {
//				int next = MathUtils.nextInt(totalElements);
//				indicator[next]=true;
//			}
			
		}
		long time2 = System.currentTimeMillis();
		System.out.println("AlignmentMapping with boolean[]:\t"+(time2 - time1) +"\t"+ (time2 - time1)/((double)(ite)) );
		
		
	}


	private static void useHashMapMapping(int ite, int swapLength, AlignmentMapping aMap) {

//		System.out.print("For count: "+totalCount +"\tElement:"+ totalElements +"\t");
		
		long time1 = System.currentTimeMillis();
		
		int[] twoPositions = new int[2];
		Set<Integer> allSrpPos = new HashSet<Integer>();
		int hapLength = aMap.getLength();
		for (int t = 0; t < ite; t++) {
			
			twoPositions[0] = MathUtils.nextInt(hapLength - swapLength);
			twoPositions[1]= twoPositions[0] + swapLength;
			
			allSrpPos.clear();
			for (int i = twoPositions[0]; i < twoPositions[1]; i++) {
				ArrayList<Integer> mapToSrp = aMap.getMapToSrp(i);
				allSrpPos.addAll(mapToSrp);
			}
			
		}
		long time2 = System.currentTimeMillis();
		System.out.println("AlignmentMapping with HashMap:\t"+(time2 - time1) +"\t"+ (time2 - time1)/((double)(ite)) );
		
		
	}

	private static AlignmentMapping setup() throws Exception{
			String dataDir = "/home/sw167/workspaceSrp/ABI/unittest/";
			String shortReadFile = "benchmark_5000Srp.fasta";
			
			DataImporter dataImporter = new DataImporter(dataDir);
			Alignment shortReads = dataImporter.importShortReads(shortReadFile);
			AlignmentMapping aMap = new AlignmentMapping(shortReads);
			return aMap;
		}

	private static void basicTest(){
		
		int ite = 10000;
//		useBoolean(ite, 100,	100);
//		useBoolean(ite, 1000,	1000);
		useBoolean(ite, 1000,	10000);
		useBoolean(ite, 1000,	100000);
//		useBoolean(ite, 10000,	10000);
//		useBoolean(ite, 10000,	100000);
//		useBoolean(ite, 10000,	1000000);
//		useBoolean(ite, 10000,	10000000);

//		useBoolean(ite, 100000,	1000);
//		useBoolean(ite, 100000,	10000);
//		useBoolean(ite, 100000,	100000);
//		
//		
//		useHashSet(ite, 100,	100);
//		useHashSet(ite, 1000,	1000);
//		useHashSet(ite, 1000,	10000);
//		useHashSet(ite, 1000,	100000);
//		useHashSet(ite, 10000,	1000000);
//		useHashSet(ite, 100000,	10000);
//		
//		reuseHashSet(ite, 100,		100);
//		reuseHashSet(ite, 1000,		1000);
//		reuseHashSet(ite, 10000,	10000);
//		reuseHashSet(ite, 10000,	100000);
//		reuseHashSet(ite, 10000,	1000000);
//		reuseHashSet(ite, 10000,	10000000);
//		reuseHashSet(ite, 10000,	1000000);
		
//		useBoolean(ite, 100000,	1000);
//		useBoolean(ite, 100000,	10000);
//		useBoolean(ite, 100000,	100000);
//		reuseHashSet(ite, 100000,	1000);
//		reuseHashSet(ite, 100000,	10000);
//		reuseHashSet(ite, 100000,	100000);
		
		
		
	}
	
	private static void useHashSet(int ite, int totalCount, int totalElements) {
		System.out.print("For count: "+totalCount +"\tElement:"+ totalElements +"\t");

		long time1 = System.currentTimeMillis();
		

		for (int t = 0; t < ite; t++) {
			HashSet<Integer> generated = new HashSet<Integer>();
			for (int i = 0; i < totalCount; i++) {
				Integer next = MathUtils.nextInt(totalElements);
				generated.add(next);
			}
			int[] uniqueArray = Ints.toArray(generated);
		}
		long time2 = System.currentTimeMillis();
		System.out.println("Use HashSet:\t"+(time2 - time1) +"\t"+ (time2 - time1)/((double)(ite)) );
		
	}


	private static void reuseHashSet(int ite, int totalCount, int totalElements) {
		
		System.out.print("For count: "+totalCount +"\tElement:"+ totalElements +"\t");
		long time1 = System.currentTimeMillis();
		HashSet<Integer> generated = new HashSet<Integer>();
		for (int t = 0; t < ite; t++) {
			generated.clear();
			for (int i = 0; i < totalCount; i++) {
				Integer next = MathUtils.nextInt(totalElements);
				generated.add(next);
			}
			int[] uniqueArray = Ints.toArray(generated);
		}
		long time2 = System.currentTimeMillis();
		System.out.println("Reuse HashSet:\t"+(time2 - time1) +"\t"+ (time2 - time1)/((double)(ite)) );

	}

	private static void useBoolean(int ite, int totalCount, int totalElements) {
		
		System.out.print("For count: "+totalCount +"\tElement:"+ totalElements +"\t");
		
		long time1 = System.currentTimeMillis();
		for (int t = 0; t < ite; t++) {
		
	
			boolean[] indicator = new boolean[totalElements];
			Arrays.fill(indicator, false);
			
			for (int i = 0; i < totalCount; i++) {
				int next = MathUtils.nextInt(totalElements);
				indicator[next]=true;
			}
			
		}
		long time2 = System.currentTimeMillis();
		System.out.println("Use boolean[]:\t"+(time2 - time1) +"\t"+ (time2 - time1)/((double)(ite)) );
		
	}
}
