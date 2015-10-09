package srp.evolution.shortreads;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;

import srp.operator.haplotypes.AbstractHaplotypeOperator;
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


/*

example mapToSrp[i].size
49	49	53	53	53	54	55	56	57	60	61	61	61	61	61	61	61	62	62	62	63	63	63	63	63	63	64	64	64	64	64	64	66	66	67	67	68	68	69	70	70	70	71	74	74	75	75	76	76	76	77	78	79	79	80	81	81	81	83	83	83	84	85	87	87	88	88	89	89	89	90	91	92	92	92	92	92	93	93	93	95	95	95	96	96	96	97	97	98	98	98	101	102	103	103	105	105	105	105	105	106	107	108	108	108	108	109	110	110	111	112	112	113	113	113	114	115	115	115	115	
116	116	116	116	116	116	116	116	117	117	117	117	117	117	117	117	118	119	119	119	120	120	120	122	123	124	124	126	126	126	127	127	128	127	127	127	128	127	126	126	126	127	127	127	127	127	127	128	127	126	127	127	127	127	128	128	128	128	128	129	130	131	129	128	130	130	131	130	130	131	132	134	135	135	133	133	133	134	135	136	136	135	134	134	134	133	132	131	131	132	130	130	130	132	131	130	130	130	130	128	128	129	128	127	130	130	130	130	129	131	132	132	132	133	132	132	134	135	134	134	
133	132	133	133	133	132	130	130	131	131	131	131	132	133	133	133	133	133	134	133	132	134	133	133	132	131	133	134	134	135	136	136	135	136	138	138	139	139	140	140	141	140	141	141	141	140	141	141	141	139	139	141	144	143	143	143	143	144	144	144	142	142	141	140	140	142	142	142	142	142	142	142	143	143	143	143	144	144	145	144	145	145	145	145	143	143	142	145	145	146	144	144	143	141	141	141	141	140	140	141	140	140	140	142	141	142	141	141	141	141	142	143	143	143	143	143	142	141	140	140	
141	142	140	140	140	140	140	140	139	139	138	140	140	140	141	141	141	142	144	145	144	143	143	143	144	144	145	146	145	143	143	143	144	145	145	145	144	145	148	148	148	147	146	145	144	144	144	144	146	149	151	151	149	150	150	153	152	152	152	154	154	154	154	156	157	156	156	158	158	159	158	156	155	156	157	156	156	156	155	154	154	154	153	152	153	151	152	154	153	153	154	154	154	153	152	151	153	153	151	150	151	151	153	155	155	153	154	153	152	151	151	151	151	150	150	149	151	152	152	152	
153	150	149	148	148	147	147	146	148	149	150	149	148	147	147	147	146	145	144	144	145	145	146	146	146	146	146	147	147	146	146	147	147	148	148	149	149	149	149	149	146	147	148	147	147	146	146	146	149	149	148	148	148	148	149	150	151	153	151	152	152	153	153	154	154	153	153	153	151	151	151	151	149	150	148	149	148	148	149	150	150	149	149	147	147	147	147	146	145	145	144	144	144	145	145	145	145	143	142	144	146	147	148	149	149	150	150	150	149	148	147	147	147	147	147	147	146	145	146	147	
147	146	148	150	150	150	150	150	150	150	150	151	151	151	151	152	151	152	153	153	154	153	152	153	154	155	155	155	156	155	154	155	155	154	154	156	155	154	153	153	154	154	155	156	156	156	157	157	155	155	155	154	154	154	153	154	153	154	155	156	154	156	154	149	149	149	148	148	149	149	148	146	147	147	148	148	148	147	147	146	144	143	143	140	140	141	141	140	139	140	142	141	141	141	141	140	140	140	140	141	141	141	141	141	141	141	142	142	143	144	144	145	145	144	145	142	143	141	142	141	
141	141	142	142	142	143	145	147	146	146	146	147	144	145	145	140	140	139	139	138	138	138	139	138	138	138	139	140	138	138	140	143	143	146	146	146	145	145	145	145	144	144	145	146	145	145	145	144	144	144	144	144	147	147	148	145	143	141	142	142	142	143	143	143	144	143	143	143	142	143	142	142	140	139	138	138	139	139	139	139	140	140	140	139	139	137	137	137	137	136	135	137	136	138	139	139	140	138	139	139	140	140	140	139	139	140	140	141	140	140	140	140	141	141	140	139	140	140	141	141	
140	140	140	141	141	141	141	142	142	140	139	140	140	141	141	141	141	140	139	138	138	138	137	140	141	140	140	140	140	138	137	137	136	136	139	139	136	138	137	136	135	136	136	136	135	136	136	137	137	137	139	139	140	140	142	142	142	142	142	142	143	143	141	141	141	140	141	140	139	139	139	139	140	140	141	142	142	142	143	145	144	144	143	143	142	144	145	144	145	146	145	146	146	146	145	143	143	142	144	142	141	142	142	141	141	141	142	142	143	144	146	147	146	145	142	142	142	143	142	142	
143	144	145	145	146	145	145	144	144	143	145	145	141	143	143	143	144	145	145	146	146	147	145	144	145	145	145	144	144	145	145	144	143	144	144	143	145	146	145	143	143	143	143	143	142	142	140	140	141	141	142	142	140	141	143	143	143	145	146	144	142	142	141	140	140	140	140	141	142	142	142	141	138	138	139	139	138	138	138	137	137	137	137	137	137	137	137	137	136	136	134	133	132	132	129	129	128	127	127	125	125	124	124	123	123	121	121	119	119	119	118	117	116	116	116	115	115	114	114	114	
113	113	113	112	112	112	112	110	110	109	109	106	104	102	102	102	99	98	98	98	98	98	97	96	96	96	96	96	94	94	94	93	90	90	90	88	88	88	88	88	88	88	87	86	86	85	84	83	83	83	83	83	83	82	81	81	80	79	79	79	77	77	77	76	75	73	73	72	72	72	72	72	72	70	69	69	68	65	65	64	64	63	63	63	62	62	62	62	61	61	61	59	59	59	59	59	59	59	58	57	54	54	54	51	50	50	50	48	46	46	45	43	43	42	41	41	40	36	23	remove taxa_0	0

every position ~ 100ish reads.
 */
public class ShortReadMapping {

	private static final int GAP = '-';
	private static final DataType DATA_TYPE = Nucleotides.INSTANCE;
	public static final char[] DNA_CHARS = AbstractHaplotypeOperator.DNA_CHARS;
	private static final double[] EQUAL_FREQ = new double[]{0.25, 0.5, 0.75, 1};
	
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
	
	private boolean isMinProp = true;
	private double minPercentage = 1; 
	
	
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

        srpCumFreqArray = new double[fullHaplotypeLength]['U'];
        
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
		createCumFreqArray();
		//TODO removed these cause old unit test to fail

	}
	
	private void createCumFreqArray() {

//		String[] srpArray = getSrpArray();
//		srpCumFreqArray = new double[fullHaplotypeLength][4];
//		for
		
		for (int i = 0; i < srpCumFreqArray.length; i++) {
			srpCumFreqArray[i]  = new double[] {
					srpCumFreqArray[i]['A'], srpCumFreqArray[i]['C'], 
					srpCumFreqArray[i]['G'], srpCumFreqArray[i]['T'] };
			
			if(isMinProp){
				
				int numZero = 0;
				double sum = 0;
				for (int j = 0; j < srpCumFreqArray[i].length; j++) {
					sum += srpCumFreqArray[i][j];
					if( srpCumFreqArray[i][j] == 0 ){
						numZero++;
					}
				}
				
				double zeroModifier = minPercentage *sum / (100-numZero*minPercentage ); 
				sum += (zeroModifier* numZero);
				for (int j = 0; j < srpCumFreqArray[i].length; j++) {
					if( srpCumFreqArray[i][j] == 0 ){
						srpCumFreqArray[i][j] = zeroModifier;
					}
					srpCumFreqArray[i][j] /= sum;
				}
				for (int j = 1; j < srpCumFreqArray[i].length; j++) {
					srpCumFreqArray[i][j] = srpCumFreqArray[i][j] + srpCumFreqArray[i][j-1];
				}
//				srpCumFreqArray[i][3]=1;
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
			System.out.println(Arrays.toString(srpCumFreqArray[i]));

		}
		
	
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
					srpCumFreqArray[j][c]++;
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
//		
		double d = MathUtils.nextDouble();
		for (int i = 0; i < 3; i++) {
			if (d <= srpCumFreqArray[pos][i]) {
				return DNA_CHARS[i];
//				return DATA_TYPE.getChar(i); //TODO: test
			}
		}
		return DNA_CHARS[3];

		//		return newChar;
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
//	public int nextBaseEqualFreq(){
//		return nextDNAFromCumFreq(EQUAL_FREQ);
//	}
//	
//	private static int nextDNAFromCumFreq(double[] cumFreq){
//		
//		double d = MathUtils.nextDouble();
//		for (int i = 0; i < cumFreq.length; i++) {
//			if (d <= cumFreq[i]) {
//				return DATA_TYPE.getChar(i); //TODO: test
//			}
//		}
//		return GAP;
//	}
//	

	
}
