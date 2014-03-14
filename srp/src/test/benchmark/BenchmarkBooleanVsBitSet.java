package test.benchmark;

import java.util.ArrayList;
import java.util.Arrays;
//import java.util.BitSet;
import java.util.HashSet;
import java.util.Set;

import srp.core.DataImporter;
import srp.shortreads.ShortReadMapping;
import cern.colt.bitvector.BitVector;

import com.carrotsearch.hppc.BitSet;
import com.google.common.primitives.Ints;

import dr.evolution.alignment.Alignment;
import dr.math.MathUtils;
/*
700 reads, ~130-150x in the middle
AlignmentMapping with boolean[]:	3274	0.003274
AlignmentMapping with HashMap:		30020	0.03002

5000 reads
AlignmentMapping with boolean[]:	3540	0.0354
AlignmentMapping with HashMap:	26402	0.26402


AlignmentMapping with BitSet:	5378	0.005378	1.958213188693E12
AlignmentMapping with BitVector:	16364	0.016364	9.56970202E8
1M reads??
*/


public class BenchmarkBooleanVsBitSet {

	public static void main(String[] args) throws Exception {

		int ite = (int) 1e6;
//		basicTest();
//		
		ShortReadMapping sMap = setup();
//		useBooleanMapping(ite, 12, sMap);
		useBitSetMapping(ite, 12, sMap);
		useBitVectorMapping(ite, 12, sMap);
//		useHashMapMapping(ite, 12, sMap);
	}
	
	private static void useBooleanMapping(int ite, int swapLength, ShortReadMapping aMap) {

//		System.out.print("For count: "+totalCount +"\tElement:"+ totalElements +"\t");
		
		long time1 = System.currentTimeMillis();
		double sum = 0;
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
			for (int i = 0; i < srpSwitch.length; i++) {
				if(srpSwitch[i]){
					sum+=i;
				}
			}
//			for (int i = 0; i < totalCount; i++) {
//				int next = MathUtils.nextInt(totalElements);
//				indicator[next]=true;
//			}
			
		}
		long time2 = System.currentTimeMillis();
		System.out.println("AlignmentMapping with boolean[]:\t"+(time2 - time1) +"\t"+ (time2 - time1)/((double)(ite)) +"\t"+ sum);
		
		
	}

	private static void useBitSetMapping(int ite, int swapLength, ShortReadMapping aMap) {

//		System.out.print("For count: "+totalCount +"\tElement:"+ totalElements +"\t");
		
		long time1 = System.currentTimeMillis();
		BitSet bs = new BitSet(aMap.getSrpCount());
		double sum = 0;
		int[] twoPositions = new int[2];
		int hapLength = aMap.getLength();
		for (int t = 0; t < ite; t++) {
			twoPositions[0] = MathUtils.nextInt(hapLength - swapLength);
			twoPositions[1]= twoPositions[0] + swapLength;

			bs.clear();
			for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
				BitSet bitSet = aMap.getBitSet(k);
				bs.or(bitSet);
			}
			int index = -1;
			sum += bs.cardinality();
			do{
				
				index = bs.nextSetBit(index+1);
				sum += index;
//				System.out.println(index);
			}while(index!= -1);
			
		}
		long time2 = System.currentTimeMillis();
		System.out.println("AlignmentMapping with BitSet:\t"+(time2 - time1) +"\t"+ (time2 - time1)/((double)(ite)) +"\t"+ sum);
		
		
	}


	private static void useBitVectorMapping(int ite, int swapLength, ShortReadMapping aMap) {

//		System.out.print("For count: "+totalCount +"\tElement:"+ totalElements +"\t");
		
		long time1 = System.currentTimeMillis();
		int srpCount = aMap.getSrpCount();
		BitVector bs = new BitVector(srpCount);
		double sum = 0;
		int[] twoPositions = new int[2];
		int hapLength = aMap.getLength();
		for (int t = 0; t < ite; t++) {
//			System.out.println(t);
			twoPositions[0] = MathUtils.nextInt(hapLength - swapLength);
			twoPositions[1]= twoPositions[0] + swapLength;
			bs.clear();
//			Arrays.fill(srpSwitch, false);
			for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
				BitVector bitSet = aMap.getBitVector(k);
				bs.or(bitSet);
				
			}
//			System.out.println(bs.toString());
			sum += bs.cardinality();
			int index2 = -1;
			do{
				index2 +=1;
				index2 = bs.indexOfFromTo(index2, 10, true);
//				System.out.println(index2);
				sum += index2;
//				index +=1;
			}while(index2 >=0 && index2<10);
			
//			do{
//				index = bs.indexOfFromTo(index+1, hapLength, true);
//				
//				sum+= index;
//			}while(index<= hapLength);
			
		}
		long time2 = System.currentTimeMillis();
		System.out.println("AlignmentMapping with BitVector:\t"+(time2 - time1) +"\t"+ (time2 - time1)/((double)(ite)) +"\t"+ sum);
		
		
	}

	private static void useHashMapMapping(int ite, int swapLength, ShortReadMapping aMap) {

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

	private static ShortReadMapping setup() throws Exception{
			String dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/";
			String shortReadFile = "benchmark_5000Srp.fasta";
			
			DataImporter dataImporter = new DataImporter(dataDir);
			Alignment shortReads = dataImporter.importShortReads(shortReadFile);
			ShortReadMapping aMap = new ShortReadMapping(shortReads);
			return aMap;
		}

	private static void basicTest(){
		
		int ite = 10000;
		useBoolean(ite, 100,	100);
		useBoolean(ite, 1000,	1000);
		useBoolean(ite, 1000,	10000);
		useBoolean(ite, 1000,	100000);
		useBoolean(ite, 10000,	10000);
		useBoolean(ite, 10000,	100000);
		useBoolean(ite, 10000,	1000000);
//		useBoolean(ite, 10000,	10000000);

//		useBoolean(ite, 100000,	1000);
//		useBoolean(ite, 100000,	10000);
//		useBoolean(ite, 100000,	100000);
//		
//		
		useBitSet(ite, 100,	100);
		useBitSet(ite, 1000,	1000);
		useBitSet(ite, 1000,	10000);
		useBitSet(ite, 1000,	100000);
		useBitSet(ite, 10000,	10000);
		useBitSet(ite, 10000,	100000);
		useBitSet(ite, 10000,	1000000);
//		useBitSet(ite, 10000,	10000000);

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
	
	private static void useBitSet(int ite, int totalCount, int totalElements) {
		System.out.print("For count: "+totalCount +"\tElement:"+ totalElements +"\t");

		long time1 = System.currentTimeMillis();
		BitSet bs = new BitSet(totalElements);
		double sum = 0;
		for (int t = 0; t < ite; t++) {
			bs.clear();
//			HashSet<Integer> generated = new HashSet<Integer>();
			for (int i = 0; i < totalCount; i++) {
				int next = MathUtils.nextInt(totalElements);
				bs.set(next);
			}
//			int[] uniqueArray = Ints.toArray(generated);
			int index = 0;
			do{
				index = bs.nextSetBit(index+1);
				sum+= index;
			}while(index!= -1);
		}
		long time2 = System.currentTimeMillis();
		System.out.println("Use BitSet:\t"+(time2 - time1) +"\t"+ (time2 - time1)/((double)(ite)) +"\t"+ sum);
		
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
		double sum = 0;
		long time1 = System.currentTimeMillis();
		for (int t = 0; t < ite; t++) {
		
			boolean[] indicator = new boolean[totalElements];
			Arrays.fill(indicator, false);
			
			for (int i = 0; i < totalCount; i++) {
				int next = MathUtils.nextInt(totalElements);
				indicator[next]=true;
			}
			for (int i = 0; i < indicator.length; i++) {
				if(indicator[i]){
					sum+=i;
				}
			}
			
		}
		long time2 = System.currentTimeMillis();
		System.out.println("Use boolean[]:\t"+(time2 - time1) +"\t"+ (time2 - time1)/((double)(ite)) +"\t"+ sum);
		
	}
}
