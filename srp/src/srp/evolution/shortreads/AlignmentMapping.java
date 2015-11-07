package srp.evolution.shortreads;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;


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

@Deprecated
public class AlignmentMapping {

	private static final int GAP = '-';
	private static final DataType DATA_TYPE =Nucleotides.INSTANCE;

	private static final double[] EQUAL_FREQ = new double[]{0.25, 0.5, 0.75, 1};
	
	private ArrayList<Integer>[] mapToSrp; // each [] = position, each ArrayList
											// = map to which read
	private HashMap<String, Integer> seqNameToSeqID; // map sequence_name >xxxto int

	private ArrayList<Character>[] listOfAvailableChar;
	private ArrayList<int[]> listOfAvailableChar2 = new ArrayList<int[]>();
	private ArrayList<ShortRead> shortReads;

	private int haplotypeLength;
	private Integer srpCount;
	private double[] cumFreq;
	

	private int[] posChar = new int[2];
	private int[] consensus;
	
	private void init(int l) {
		haplotypeLength = l;
		srpCount = 0;

		cumFreq = new double['U'];

		seqNameToSeqID = new HashMap<String, Integer>();
		shortReads = new ArrayList<ShortRead>();

		listOfAvailableChar = new ArrayList[haplotypeLength];
		mapToSrp = new ArrayList[haplotypeLength];
		for (int i = 0; i < haplotypeLength; i++) {
			mapToSrp[i] = new ArrayList<Integer>();
		}

	}
	@Deprecated
	public AlignmentMapping(Alignment srpAlignment){
		
		init( srpAlignment.getSiteCount() );
		
		@SuppressWarnings("unchecked")
		HashSet<Character>[] setsOfAvailableChar = new HashSet[haplotypeLength];

		for (int i = 0; i < setsOfAvailableChar.length; i++) {
			setsOfAvailableChar[i] = new HashSet<Character>();
		}
		
        int[][] frequencies = new int[haplotypeLength][DATA_TYPE.getAmbiguousStateCount()];


		for (int i = 0; i < srpAlignment.getSequenceCount(); i++) {
			Sequence s = srpAlignment.getSequence(i);
			addSequence(s, setsOfAvailableChar, frequencies);
		}
		//post processing after created the basic framework
		createConsensusSequence(frequencies);
    
		//calculated listOfAvailableChar(2)
		calculateListOfAvailableChar(setsOfAvailableChar);
		calculateCumFreq();
		//TODO removed these cause old unit test to fail

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

	@Deprecated
	private void calculateListOfAvailableChar(HashSet<Character>[] setsOfAvailableChar) {
		for (int i = 0; i < haplotypeLength; i++) {
			listOfAvailableChar[i] = new ArrayList<Character>(setsOfAvailableChar[i]);

			Character[] temp = new Character[setsOfAvailableChar[i].size()];
			setsOfAvailableChar[i].toArray(temp);
			int[] temp2 = new int[temp.length];
			for (int j = 0; j < temp2.length; j++) {
				temp2[j] = temp[j];
			}
			
			listOfAvailableChar2.add(temp2); //TODO: even faster??
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
		return haplotypeLength;
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
		return shortReads.get(i).getLength();
	}

	public String getSrpName(int i) {
		return shortReads.get(i).getName();
	}


	protected int nextBaseAt(int pos){
		int newChar = GAP;
		int size = mapToSrp[pos].size();
		if (size != 0) {
			int srpIndex = mapToSrp[pos].get(MathUtils.nextInt(size));
			newChar = getShortReadCharAt(srpIndex, pos);
		}
		
		return newChar;
	}


	public int[] getNextBase() {

		int pos = MathUtils.nextInt(haplotypeLength);
		int newChar = nextBaseAt(pos);
	
		return new int[]{pos, newChar};
	}
	

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
	
	 
	  
	 
	public int[] getNextBaseUniform() {
		posChar[0] = MathUtils.nextInt(haplotypeLength);
		int[] chars = listOfAvailableChar2.get(  posChar[0] );
		int size = chars.length;
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
		if (size != 0) {
			posChar[1] = chars[ MathUtils.nextInt(size) ];
		}
		else{
			posChar[1] = GAP;
		}
		
		return posChar;
	}
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
	public int[] getNextBaseEmpirical() {

		posChar[0] = MathUtils.nextInt(haplotypeLength);
		
//		double d = MathUtils.nextDouble();
//
//		for (int i = 0; i < cumFreq.length; i++) {
//			if (d <= cumFreq[i]) {
//				posChar[1] = VALID_CHARS[i];
//				return posChar;
//			}
//		}
//		posChar[1] = GAP;
		posChar[1] = nextDNAFromCumFreq(cumFreq);
		return posChar;
	}
	
	public int nextBaseEqualFreq(){
		return nextDNAFromCumFreq(EQUAL_FREQ);
	}
	
	private static int nextDNAFromCumFreq(double[] cumFreq){
		
		double d = MathUtils.nextDouble();
		for (int i = 0; i < cumFreq.length; i++) {
			if (d <= cumFreq[i]) {
				return DATA_TYPE.getChar(i); //TODO: test
			}
		}
		return GAP;
	}
	

	
}
