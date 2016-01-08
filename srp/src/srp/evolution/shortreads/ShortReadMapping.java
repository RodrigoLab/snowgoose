package srp.evolution.shortreads;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.lang3.BooleanUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.ArithmeticUtils;

import srp.operator.haplotypes.AbstractHaplotypeOperator;
import srp.distributions.DirichletMultinomialDistribution;
import srp.likelihood.AbstractShortReadsLikelihood;
import cern.colt.bitvector.BitVector;

import com.google.common.primitives.Ints;
import com.sun.org.apache.bcel.internal.generic.DNEG;

import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;
//import java.util.BitSet;


public class ShortReadMapping {
	
	private static final int GAP = '-';
	private static final DataType DATA_TYPE = Nucleotides.INSTANCE;
	private static final double[] EQUAL_FREQ = new double[]{0.25, 0.5, 0.75, 1};
	
	public static final char[] DNA_CHARS = AbstractHaplotypeOperator.DNA_CHARS;
	public static final double LOG_ERROR_RATE = AbstractShortReadsLikelihood.LOG_ERROR_RATE;
	public static final double LOG_ONE_MINUS_ERROR_RATE = AbstractShortReadsLikelihood.LOG_ONE_MINUS_ERROR_RATE;

	
	
	private static final double MIN_PERCENTAGE = 1;//0.1; 
	private static final int MIN_READ_COUNT = 2;
	private static final boolean IS_MIN_PROPORTION = true;
	private static final double LOW_VARIANCE_THRESHOLD = 0.05;
	
	private ArrayList<Integer>[] mapToSrp; // each [] = position, each ArrayList map to which read
	
	private HashMap<String, Integer> seqNameToSeqID; // map sequence_name >xxxto int

	private ArrayList<Character>[] listOfAvailableChar;
	private ArrayList<int[]> listOfAvailableChar2 = new ArrayList<int[]>();
	private ArrayList<ShortRead> shortReads;

	private int fullHaplotypeLength;
	private Integer srpCount;
	private double[] cumFreq;
	

//	private int[] posChar = new int[2];
	private int[] consensus;
	private String[] srpArray;
	private int[][] mapToSrpArray;
	private BitSet[] bitSetArray;
	private BitVector[] bitVectorArray;
	private int[] srpLength;
	private char[][] srpChar2D;
	private Integer[] allSrpLengthInteger;
	
	private double[][] srpCumFreqArray;
	private int[][] srpCountArray;
	

	private boolean[] isFixedSite;
	private boolean[] isLowVarianceSite;
	private int lvs;
	private SiteType[] siteTypesArray;
	
	private void init(int l) {
		fullHaplotypeLength = l;
		srpCount = 0;

		cumFreq = new double['U'];

		seqNameToSeqID = new HashMap<String, Integer>();
		shortReads = new ArrayList<ShortRead>();

		listOfAvailableChar = new ArrayList[fullHaplotypeLength];
		mapToSrp = new ArrayList[fullHaplotypeLength];
		for (int i = 0; i < fullHaplotypeLength; i++) {
			mapToSrp[i] = new ArrayList<Integer>();
		}

	}

	public ShortReadMapping(Alignment srpAlignment){
		
		init( srpAlignment.getSiteCount() );
		
		@SuppressWarnings("unchecked")
		HashSet<Character>[] setsOfAvailableChar = new HashSet[fullHaplotypeLength];

		for (int i = 0; i < setsOfAvailableChar.length; i++) {
			setsOfAvailableChar[i] = new HashSet<Character>();
		}
		
        int[][] frequencies = new int[fullHaplotypeLength][DATA_TYPE.getAmbiguousStateCount()];

        srpCountArray = new int[fullHaplotypeLength]['T'+1];
        srpCumFreqArray = new double[fullHaplotypeLength][4];
        
		for (int i = 0; i < srpAlignment.getSequenceCount(); i++) {
			Sequence s = srpAlignment.getSequence(i);
			addSequence(s, setsOfAvailableChar, frequencies);
		}
		//post processing after created the basic framework
		createConsensusSequence(frequencies);

		createSrpArray();
		createMapToSrpArray();
		createBitSetArray();
		createSrpChar2DArray();
		//calculated listOfAvailableChar(2)
		createListOfAvailableChar(setsOfAvailableChar);
		calculateCumFreq();
		createCumFreqArray(); //Use counts to create prob
//		createProbForEachBase(); // likelihood -> prob, mdm
		//TODO removed these cause old unit test to fail
		System.out.println(lvs);

	}
	private void createProbForEachBase(){
		DirichletMultinomialDistribution dmd = new DirichletMultinomialDistribution();
		
		for (int i = 0; i < srpCumFreqArray.length; i++) {
			
//			srpCumFreqArray[i]['A'], srpCumFreqArray[i]['C'], 
//			srpCumFreqArray[i]['G'], srpCumFreqArray[i]['T'] 
			double[] probBases = new double[4];
			int sum = 	srpCountArray[i]['A'] + srpCountArray[i]['C'] + 
					srpCountArray[i]['G'] + srpCountArray[i]['T'] ;
			
			probBases[0] = ArithmeticUtils.binomialCoefficientLog(sum, srpCountArray[i]['A'])+
					srpCountArray[i]['A']*LOG_ONE_MINUS_ERROR_RATE+(sum-srpCountArray[i]['A'])*LOG_ERROR_RATE;
			
			probBases[1]  = ArithmeticUtils.binomialCoefficientLog(sum, srpCountArray[i]['C'])+
					srpCountArray[i]['C']*LOG_ONE_MINUS_ERROR_RATE+(sum-srpCountArray[i]['C'])*LOG_ERROR_RATE;
					
			probBases[2] = ArithmeticUtils.binomialCoefficientLog(sum, srpCountArray[i]['G'])+
					srpCountArray[i]['G']*LOG_ONE_MINUS_ERROR_RATE+(sum-srpCountArray[i]['G'])*LOG_ERROR_RATE;
					
			probBases[3] = ArithmeticUtils.binomialCoefficientLog(sum, srpCountArray[i]['T'])+
					srpCountArray[i]['T']*LOG_ONE_MINUS_ERROR_RATE+(sum-srpCountArray[i]['T'])*LOG_ERROR_RATE;

//			System.out.println(srpCountArray[i]['A'] +"\t"+ srpCountArray[i]['C'] +"\t"+ srpCountArray[i]['G'] +"\t"+ srpCountArray[i]['T']);
			System.out.println(Arrays.toString(probBases));
			
			int[] counts =  new int[] {
					srpCountArray[i]['A'], srpCountArray[i]['C'], 
					srpCountArray[i]['G'], srpCountArray[i]['T'] };
//			counts =  new int[] {
//					0,50,49,51};

			double[] result = dmd.FourDirichletMultinomialLogProbability(counts);
//			double max = StatUtils.max(result);
//			for (int j = 0; j < result.length; j++) {
//				probBases[j] = result[j] - max;
//			}
			probBases = result;
			System.out.println(Arrays.toString(probBases));
			double sumProb = 0;
			for (int j = 0; j < probBases.length; j++) {
				probBases[j] = Math.exp(probBases[j]);
				sumProb += probBases[j];
			}
//			System.out.println(sumProb);
//			
			for (int j = 0; j < probBases.length; j++) {
				probBases[j] /= sumProb;
			}
//			System.out.println(Arrays.toString(probBases));
			srpCumFreqArray[i][0] = probBases[0];
			for (int j = 1; j < probBases.length; j++) {
				srpCumFreqArray[i][j] = probBases[j] + srpCumFreqArray[i][j-1];
			}
			if(srpCumFreqArray[i][3]!= 1){
				srpCumFreqArray[i][3] = 1;
				
			}
			
			System.out.println(Arrays.toString(result));
			System.out.println(srpCountArray[i]['A'] +"\t"+ srpCountArray[i]['C'] +"\t"+ srpCountArray[i]['G'] +"\t"+ srpCountArray[i]['T']);
			System.out.println(Arrays.toString(probBases));
//			System.out.println("srpCumFreqArray: "+Arrays.toString(srpCumFreqArray[i]));
			System.out.println();
		}

		
		System.exit(12);
		
	}
	private void createCumFreqArray() {

//		String[] srpArray = getSrpArray();
//		srpCumFreqArray = new double[fullHaplotypeLength][4];
//		for
//		DirichletMultinomialDistribution dmd = new DirichletMultinomialDistribution();
		isFixedSite = new boolean[srpCumFreqArray.length];
		isLowVarianceSite = new boolean[srpCumFreqArray.length];
		siteTypesArray = new SiteType[srpCumFreqArray.length];
		
		for (int i = 0; i < srpCumFreqArray.length; i++) {
			srpCumFreqArray[i]  = new double[] {
					srpCountArray[i]['A'], srpCountArray[i]['C'], 
					srpCountArray[i]['G'], srpCountArray[i]['T'] };
//			System.out.println(Arrays.toString(srpCumFreqArray[i]));
			if(IS_MIN_PROPORTION){
				
				int numMinCount = 0;
				double sum = 0;
				for (int j = 0; j < srpCumFreqArray[i].length; j++) {
					sum += srpCumFreqArray[i][j];
					if( srpCumFreqArray[i][j] <= MIN_READ_COUNT ){
						numMinCount++;
					}
				}

				
				double minModifier = MIN_PERCENTAGE *sum / (100-numMinCount*MIN_PERCENTAGE );
				
				if(numMinCount == 3){
					isFixedSite[i] = true;
					siteTypesArray[i] = SiteType.FIXED;
//					minModifier = 0;
				}
//				System.out.println(isFixedSite[i]);
//				sum += (zeroModifier* numMinCount);
				sum = 0;
				for (int j = 0; j < srpCumFreqArray[i].length; j++) {
					if( srpCumFreqArray[i][j] <= MIN_READ_COUNT  ){
						srpCumFreqArray[i][j] = minModifier;
					}
					sum += srpCumFreqArray[i][j];
				}
				
				int lowVarianceCount = 0;
				if(sum == 0){
					Arrays.fill(srpCumFreqArray[i],0.25);
				}
				else{
					for (int j = 0; j < srpCumFreqArray[i].length; j++) {
						srpCumFreqArray[i][j] /= sum;
						if (srpCumFreqArray[i][j] < LOW_VARIANCE_THRESHOLD){
							lowVarianceCount++;
						}
					}
				}
				if (lowVarianceCount > 2 ){
					isLowVarianceSite[i] = true;
					siteTypesArray[i] = SiteType.LOW_VAR;
//					System.out.println(">2: "+ lowVarianceCount +"\t"+ Arrays.toString(srpCumFreqArray[i]));
					lvs++;
				}
				if (lowVarianceCount < 2 ){
					siteTypesArray[i] = SiteType.HIGH_VAR;
//					System.out.println("<2: "+ lowVarianceCount +"\t"+ Arrays.toString(srpCumFreqArray[i]));
				}
				if (lowVarianceCount == 2 ){
					siteTypesArray[i] = SiteType.HIGH_VAR;
//					System.out.println("=2: "+ lowVarianceCount +"\t"+ Arrays.toString(srpCumFreqArray[i]));
				}
				if(numMinCount == 3){
					isFixedSite[i] = true;
					siteTypesArray[i] = SiteType.FIXED;
				}
				
//				int[] counts =  new int[] {srpCountArray[i]['A'], srpCountArray[i]['C'], srpCountArray[i]['G'], srpCountArray[i]['T'] };
//				System.out.println(Arrays.toString(counts));
//				System.out.println("srpCumFreqArrayV1: "+Arrays.toString(srpCumFreqArray[i]));
//				System.out.println(isLowVarianceSite[i] +"\t"+ lowVarianceCount);
//				System.out.println();
//				
				for (int j = 1; j < srpCumFreqArray[i].length; j++) {
					srpCumFreqArray[i][j] = srpCumFreqArray[i][j] + srpCumFreqArray[i][j-1];
				}
//				srpCumFreqArray[i][3]=1;
//				System.out.println("srpCumFreqArrayV1: "+Arrays.toString(srpCumFreqArray[i]));
			}
			else{
				
				for (int j = 1; j < srpCumFreqArray[i].length; j++) {
					srpCumFreqArray[i][j] = srpCumFreqArray[i][j] + srpCumFreqArray[i][j-1];
				}
				double sum = srpCumFreqArray[i][3];
				
				for (int j = 0; j < srpCumFreqArray[i].length; j++) {
					srpCumFreqArray[i][j] /= sum;  
				}
			}
//			int[] counts =  new int[] {srpCountArray[i]['A'], srpCountArray[i]['C'], srpCountArray[i]['G'], srpCountArray[i]['T'] };
//System.out.println(Arrays.toString(counts));
//			System.out.println("srpCumFreqArrayV1: "+Arrays.toString(srpCumFreqArray[i]));
		}
//		int count1 = 0;
//		int count2 = 0;
//		for (int j = 0; j < isFixedSite.length; j++) {
//			if(isFixedSite[j]) count1++;
//			if(isLowVarianceSite[j]) count2++;
//		}
//		System.out.println("isXXSite:" + count1 +"\t"+ count2);
	
	}

	private void createBitSetArray() {
		bitSetArray = new BitSet[mapToSrp.length];
		bitVectorArray = new BitVector[mapToSrp.length];
		for (int i = 0; i < bitSetArray.length; i++) {
			bitSetArray[i] = new BitSet(srpCount);
			bitVectorArray[i] = new BitVector(srpCount);
			for (int s : mapToSrp[i]) {
				bitSetArray[i].set(s);
				bitVectorArray[i].set(s);
			}
		}
		
	}
	public BitSet getBitSet(int i){
		return bitSetArray[i];
	}
	public BitVector getBitVector(int i){
		return bitVectorArray[i];
	}
	public int[][] getSrpState2DArray() {
		String[] srpArray = getSrpArray();
		int[][] state2DArray = new int[srpArray.length][fullHaplotypeLength];
	
		for (int i = 0; i < srpArray.length; i++) {
			String srp = srpArray[i];
			for (int j = 0; j < fullHaplotypeLength; j++) {
				
				char srpChar = srp.charAt(j);//TODO: change to char[] at hight lv or at ShortReadMapping
				int state = DATA_TYPE.getState(srpChar);
				state2DArray[i][j] = state;
			}
		}
	
		return state2DArray;
	}

	public char[][] getSrpChar2DArray() {
		return srpChar2D;
	}
	
	private void createSrpChar2DArray(){
	
		String[] srpArray = getSrpArray();
		srpChar2D = new char[srpArray.length][fullHaplotypeLength];
	
		for (int i = 0; i < srpArray.length; i++) {
			String srp = srpArray[i];
			for (int j = 0; j < fullHaplotypeLength; j++) {
//				allSrpState2D[i][j] = getStateAtK(srp, j);
				srpChar2D[i][j] = srp.charAt(j);
			}
		}
		
	}

	public Integer[] getAllSrpLengthInteger() {
		return allSrpLengthInteger;
	}

	public int[][] getMapToSrpArray(){
		return mapToSrpArray;
	}

	public String[] getSrpArray(){
		return srpArray;
	}

	private void createSrpArray() {
		srpArray = new String[srpCount];
		srpLength = new int[srpCount];
		allSrpLengthInteger = new Integer[srpCount];
		for (int i = 0; i < srpArray.length; i++) {
			ShortRead shortRead = shortReads.get(i);
			srpArray[i] = shortRead.getFullSrp();
			srpLength[i] = shortRead.getLength();
			allSrpLengthInteger[i] = (Integer) shortRead.getLength();
		}
		
		
	}
	private void createMapToSrpArray(){
		mapToSrpArray = new int[mapToSrp.length][];
		for (int i = 0; i < mapToSrpArray.length; i++) {
			ArrayList<Integer> arrayList = mapToSrp[i];
			Integer[] temp = new Integer[arrayList.size()];
//			temp = arrayList.toArray(new double[arrayList.size()]);
			arrayList.toArray(temp);
			mapToSrpArray[i] = Ints.toArray(arrayList);
			
			String s1 = arrayList.toString();
			String s2 = Arrays.toString(mapToSrpArray[i]);
			if(!s1.equals(s2)){
				System.out.println(s1);
				System.out.println(s2);
			}
			
		}

		
		
		
	}
	@Deprecated
	private void calculateCumFreq() {
		cumFreq = new double[] { cumFreq['A'], cumFreq['C'], cumFreq['G'], cumFreq['T'] };
		for (int i = 1; i < cumFreq.length; i++) {
			cumFreq[i] = cumFreq[i] + cumFreq[i-1];  
		}
		double sum = cumFreq[3];
		for (int i = 0; i < cumFreq.length; i++) {
			cumFreq[i] /= sum;  
		}
		
	}

//	@Deprecated
	private void createListOfAvailableChar(HashSet<Character>[] setsOfAvailableChar) {
		for (int i = 0; i < fullHaplotypeLength; i++) {
			listOfAvailableChar[i] = new ArrayList<Character>(setsOfAvailableChar[i]);

			Character[] temp = new Character[setsOfAvailableChar[i].size()];
			setsOfAvailableChar[i].toArray(temp);
			int[] temp2 = new int[temp.length];
			for (int j = 0; j < temp2.length; j++) {
				temp2[j] = temp[j];
			}
			
			listOfAvailableChar2.add(temp2); 
//			System.out.println(listOfAvailableChar[i].toString());
//			System.out.println(Arrays.toString(  listOfAvailableChar2.get(i)  ));
//			System.out.println();
		}
	}

	private void createConsensusSequence(int[][] frequencies) {

        int[] counts;
        
        counts = new int[getLength()];
        consensus = new int[getLength()];
        for (int i = 0; i < frequencies.length; i++) {
            int maxState = 0;
            int maxFreq = frequencies[i][0];
            for (int j = 1; j < frequencies[i].length; j++) {
                int freq = frequencies[i][j];
                if (freq > maxFreq) {
                    maxState = j;
                    maxFreq = freq;
                }
            }
            consensus[i] = maxState;
            counts[i] = maxFreq;
        }

        StringBuffer buffer = new StringBuffer();
        for (int i = 0; i < consensus.length; i++) {
            buffer.append(DATA_TYPE.getChar(    consensus[i] ));
        }
        Sequence sequence = new Sequence(new Taxon("con"),buffer.toString());
        sequence.setDataType(DATA_TYPE);
//        System.out.println();
//        System.err.println(sequence.getSequenceString());
//        System.out.println();
		
	}
	private void addSequence(Sequence s,
			HashSet<Character>[] setsOfAvailableChar, 
			int[][] frequencies) {
		
		ShortRead srp = new ShortRead(s);
		if (srp.getIsValid()){
			seqNameToSeqID.put(srp.getName(), srpCount);
			shortReads.add(srpCount, srp);

			for (int j = srp.getStart(); j < srp.getEnd(); j++) {
				mapToSrp[j].add(srpCount);
				char c = srp.getFullSrpCharAt(j);
				int state = DATA_TYPE.getState(c);
				if (c!= GAP){
					setsOfAvailableChar[j].add(c);
					frequencies[j][state] += 1;
					cumFreq[c]++;
					srpCountArray[j][c]++;
				}
			}
			srpCount++;
		}
		


	}

	@Override
	public String toString() {
		String s = "";
		for (int i = 0; i < mapToSrp.length; i++) {
			s += mapToSrp[i].toString();
			s += "\n";

		}
		return s;

	}

	public int[] getConsensusSequenceState(){
		return consensus;
	}

	public ArrayList<Integer> getMapToSrp(int pos) {
		return mapToSrp[pos];
	}
	public Integer mapNameToID(String name){
		return seqNameToSeqID.get(name);
	}

	public int getLength() {
		return fullHaplotypeLength;
	}

	public int getSrpCount() {
		return srpCount;
	}

	public ShortRead getShortRead(int i) {
		return shortReads.get(i);
	}

	public char getShortReadCharAt(int index, int c) {
		return shortReads.get(index).getFullSrpCharAt(c);
	}

	public String getSrpFull(int i) {
//		return srpArray[i];//TESTING
		return shortReads.get(i).getFullSrp();
	}

	public String getSrpFragment(int i) {

		return shortReads.get(i).getFragmentSrp();
	}

	public int getSrpStart(int i) {

		return shortReads.get(i).getStart();
	}

	public int getSrpEnd(int i) {

		return shortReads.get(i).getEnd();
	}

	public int getSrpLength(int i) {
//		return shortReads.get(i).getLength();
		return srpLength[i];
	}

	public String getSrpName(int i) {
		return shortReads.get(i).getName();
	}
	public ArrayList<int[]> getListOfAvailableChar2(){
		return listOfAvailableChar2;
	}
	
	public char getBaseAt(int s) {
//		posChar[0] = MathUtils.nextInt(haplotypeLength);
		int[] chars = listOfAvailableChar2.get(  s );
		int size = chars.length;
//		System.out.println(posChar[0] +"\t"+ Arrays.toString(chars) +"\t"+ size +"\t"+ posChar[1]);
		char newChar;
//		if(size < 2){
//			return '-';
//		}
		
		if (size != 0) {
			
			newChar = (char) chars[ MathUtils.nextInt(size) ];
//			System.out.println(newChar +"\t"+ size +"\t"+ Arrays.toString(chars));
//			newChar= DNA_CHARS[ chars[ MathUtils.nextInt(size) ] ];
		}
		else{
			newChar = DNA_CHARS[ MathUtils.nextInt(4) ];
		}
		
		return newChar;
//		int[] srpList = mapToSrpArray[s];
////		System.out.println(s +"\t"+ Arrays.toString(srpList));
//		int state = 999;
//		char newChar;
//		do{
//			int srpIndex = srpList[MathUtils.nextInt(srpList.length)];
//			newChar = srpChar2D[srpIndex][s];
//			state = Nucleotides.INSTANCE.getState(newChar);
//			System.out.println(s +"\t"+ newChar +"\t"+ state);
//		}while(state>3);
//		
//		return newChar;
	}
	
	public char[] getRandHaplotype() {
		
		char[] randChar = new char[fullHaplotypeLength];
		for (int i = 0; i < randChar.length; i++) {
			randChar[i] = (char) nextBaseEqualFreq();
		}
		return randChar;
	}
	
	public char[] getSemiRandHaplotype2() {
		
		char[] randChar = new char[fullHaplotypeLength];
		for (int i = 0; i < randChar.length; i++) {
			randChar[i] = nextBaseFreqAt(i);
		}
		return randChar;
	}
	
	public char[] getSemiRandHaplotype() {
		
//		System.out.println("SemiRandHap");
		double switchSrpProb = 0.1;
		char[] randChar = new char[fullHaplotypeLength];
		
		int[] srpList = mapToSrpArray[0];
		

		int state = 999;
		char newChar ='-';
		int srpIndex = 0;
		
		if(srpList.length==0){
			newChar = 'A';//FIXME: later, only accept ACGT no N at moment
		}
		else{
			srpIndex = srpList[MathUtils.nextInt(srpList.length)];
			newChar = srpChar2D[srpIndex][0];
		}
		randChar[0] = newChar;
		
		//FIXME: not very smart loop, just get things working
		for (int s = 1; s < randChar.length; s++) {
			int count = 0;
			if(MathUtils.nextDouble() < switchSrpProb){
					srpList = mapToSrpArray[s];
					if(srpList.length!=0){
						srpIndex = srpList[MathUtils.nextInt(srpList.length)];
					}
//					state = 999;
//					System.out.println(s +"\tnewSrp Prob: "+ srpIndex);
			}
			
			do{
				if(state>3){
					if(count==10){
						newChar = 'A';
					}
					else{
						srpList = mapToSrpArray[s];
						if(srpList.length==0){
							newChar = 'A';//FIX later, only accept ACGT no N at moment
						}
						else{
							srpIndex = srpList[MathUtils.nextInt(srpList.length)];
							newChar = srpChar2D[srpIndex][s];
							count++;
						}
					}
//					System.out.println(s +"\tnewSrp: "+ srpIndex);
				}
				else{
					newChar = srpChar2D[srpIndex][s];
				}
				state = Nucleotides.INSTANCE.getState(newChar);
				
			}
			while(state>3);
			randChar[s] = newChar;
			
		}
		
		return randChar;
	}

	public void summary() {
//		System.out.println(toString());
		System.out.println(Arrays.toString(cumFreq));
//		System.out.println(Arrays.toString(consensus));
//		System.out.println(Arrays.toString(srpArray));
//		System.out.println(Arrays.toString(mapToSrpArray[0]));
		System.out.println(Arrays.toString(srpChar2D[0]));
		Integer[] x = new Integer[100];
		ArrayList<Character> aa = listOfAvailableChar[100];
		for (Character character : aa) {
			System.out.println(character);
		}
//		System.out.println(Arrays.toString(x));
		
		int[] aa2 = listOfAvailableChar2.get(100);
		System.out.println(Arrays.toString(aa2));
		
	}

//
	public int nextBaseAt(int pos){
		int newChar = GAP;
		int size = mapToSrp[pos].size();
		if (size != 0) {
			int srpIndex = mapToSrp[pos].get(MathUtils.nextInt(size));
			newChar = getShortReadCharAt(srpIndex, pos);
		}
		
		return newChar;
	}
	public char nextBaseFreqAt(int pos){
//		int newChar = GAP;
//		int size = mapToSrp[pos].size();
//		if (size != 0) {
//			int srpIndex = mapToSrp[pos].get(MathUtils.nextInt(size));
//			newChar = getShortReadCharAt(srpIndex, pos);
//		}
//		System.out.println(siteTypesArray[pos] +"\t"+ Arrays.toString(srpCumFreqArray[pos]));

		double d = MathUtils.nextDouble();
		for (int i = 0; i < 3; i++) {
			if (d <= srpCumFreqArray[pos][i]) {
				return DNA_CHARS[i];
			}
		}
		return DNA_CHARS[3];

	}

	public boolean isFixSite(int siteIndex) {
//		if(isFixedSite[siteIndex]){
//			System.out.println("isFixed: "+Arrays.toString(srpCumFreqArray[siteIndex]));
//		}
		return isFixedSite[siteIndex];
	}

	public boolean isLowVarianceSite(int siteIndex) {
//		if(isLowVarianceSite[siteIndex]){
//			System.out.println("isLowV: "+Arrays.toString(srpCumFreqArray[siteIndex]));
//		}
		return isLowVarianceSite[siteIndex];
	}

	public SiteType getSiteType(int siteIndex) {
//		System.out.println(siteTypesArray[siteIndex] +"\t"+ Arrays.toString(srpCumFreqArray[siteIndex]));
		return siteTypesArray[siteIndex];
	}
	
//
//	public int[] getNextBase() {
//
//		int pos = MathUtils.nextInt(haplotypeLength);
//		int newChar = nextBaseAt(pos);
//	
//		return new int[]{pos, newChar};
//	}
	

//	public int[] getNextBaseUniform() {
//		int newChar = GAP ;
//		int pos = MathUtils.nextInt(haplotypeLength);
//		ArrayList<Character> charList = listOfAvailableChar[pos];
//		int size = charList .size();
//		if (size != 0) {
//			newChar = charList .get(MathUtils.nextInt(size));
//		}
//		return new int[]{pos, newChar};
//	}
	
	 
	  
	 
//	public int[] getNextBaseUniform() {
//		posChar[0] = MathUtils.nextInt(haplotypeLength);
//		int[] chars = listOfAvailableChar2.get(  posChar[0] );
//		int size = chars.length;
//		System.out.println(posChar[0] +"\t"+ Arrays.toString(chars) +"\t"+ size +"\t"+ posChar[1]);
//		switch (size) {
//		case 0:
//			posChar[1] = GAP;
//			break;
//		case 1:
//			posChar[1] = chars[0];
//			break;
//		case 2:
//			boolean index = MathUtils.nextBoolean();
//			if(index){
//				posChar[1] = chars[0];
//			}
//			else{
//				posChar[1] = chars[1];
//			}
//			break;
//			
//		default:
//			posChar[1] = chars[ MathUtils.nextInt(size) ];
//			break;
//		}		
//		if (size != 0) {
//			posChar[1] = chars[ MathUtils.nextInt(size) ];
//		}
//		else{
//			posChar[1] = GAP;
//		}
//		
//		return posChar;
//	}
//	 
//	public int[] getNextBaseUniform() {
//		int pos = MathUtils.nextInt(haplotypeLength);
//		posChar[0] = pos;//MathUtils.nextInt(haplotypeLength)
//		
//		int size = listOfAvailableChar[pos].size();
//		if (size != 0) {
//			posChar[1] = listOfAvailableChar[pos].get(MathUtils.nextInt(size));
//		}
//		else{
//			posChar[1] = GAP;
//		}
//		return posChar;
//	}
//	

//	private int[] getNextBaseEqualFreq() {
//
//		posChar[0] = MathUtils.nextInt(haplotypeLength);
//		posChar[1] = nextDNAFromCumFreq(EQUAL_FREQ );
//		return posChar;
//	}
//	public int[] getNextBaseEmpirical() {
//
//		posChar[0] = MathUtils.nextInt(haplotypeLength);
//		
//		double d = MathUtils.nextDouble();
//
//		for (int i = 0; i < cumFreq.length; i++) {
//			if (d <= cumFreq[i]) {
//				posChar[1] = VALID_CHARS[i];
//				return posChar;
//			}
//		}
//		posChar[1] = GAP;
//		posChar[1] = nextDNAFromCumFreq(cumFreq);
//		return posChar;
//	}
//	
	public int nextBaseEqualFreq(){
		return nextDNAFromCumFreq(EQUAL_FREQ);
	}
//	
	private static int nextDNAFromCumFreq(double[] cumFreq){
		
		double d = MathUtils.nextDouble();
		for (int i = 0; i < cumFreq.length; i++) {
			if (d <= cumFreq[i]) {
				return DATA_TYPE.getChar(i); //TODO: test
			}
		}
		return GAP;
	}
//	

	
}
