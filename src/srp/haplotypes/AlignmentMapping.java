package srp.haplotypes;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

public class AlignmentMapping {

	private static final int GAP = '-';
	private static final DataType dataType =Nucleotides.INSTANCE;
//	private static final char[] VALID_CHARS = Nucleotides.NUCLEOTIDE_CHARS;
//	char[] x = Nucleotides.INSTANCE.getValidChars();
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

	public AlignmentMapping(Alignment srpAlignment){
		
		init( srpAlignment.getSiteCount() );
		
		@SuppressWarnings("unchecked")
		HashSet<Character>[] setsOfAvailableChar = new HashSet[haplotypeLength];

		
        int[][] frequencies = new int[haplotypeLength][dataType.getAmbiguousStateCount()];

		for (int i = 0; i < setsOfAvailableChar.length; i++) {
			setsOfAvailableChar[i] = new HashSet<Character>();
		}

		for (int i = 0; i < srpAlignment.getSequenceCount(); i++) {
			Sequence s = srpAlignment.getSequence(i);
			addSequence(s, setsOfAvailableChar, frequencies);
		}

		createConsensusSequence(frequencies);
    
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
		
		cumFreq = new double[] { cumFreq['A'], cumFreq['C'], cumFreq['G'], cumFreq['T'] };
		for (int i = 1; i < cumFreq.length; i++) {
			cumFreq[i] = cumFreq[i] + cumFreq[i-1];  
		}
		double sum = cumFreq[3];
		for (int i = 0; i < cumFreq.length; i++) {
			cumFreq[i] /= sum;  
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
            buffer.append(dataType.getChar(    consensus[i] ));
        }
        Sequence sequence = new Sequence(new Taxon("con"),buffer.toString());
        sequence.setDataType(dataType);
//        System.out.println();
//        System.err.println(sequence.getSequenceString());
//        System.out.println();
		
	}
	public int[] getConsensusSequenceState(){
		return consensus;
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
				int state = dataType.getState(c);
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
				return dataType.getChar(i); //TODO: test
			}
		}
		return GAP;
	}
	


	@Deprecated
	public int[] getNextBaseFrequency(Parameter freqs) {
//		// TODO Auto-generated method stub
//
//		posChar[0] = MathUtils.nextInt(haplotypeLength);
//		double d = MathUtils.nextDouble();
//		
////		double[] freqsValue = freqs.getParameterValues();
////		double sumFreq = 0;
////		for (int i = 0; i < freqsValue.length; i++) {
////			sumFreq += freqsValue[i];
////			if (d <= sumFreq) {
////				posChar[1] = VALID_CHARS[i];
//////				System.out.println(d +"\t"+ posChar[1] +"\t"+  sumFreq +"\t"+ Arrays.toString(freqsValue)				);
////				return posChar;
////			}
////			
////		}
////		
//		
//		double sumFreq = 0;
//		for (int i = 0; i < freqs.getDimension(); i++) {
//			sumFreq += freqs.getParameterValue(i);
//			if (d <= sumFreq) {
//				posChar[1] = VALID_CHARS[i];
////				System.out.println(d +"\t"+ posChar[1] +"\t"+  sumFreq +"\t"+ Arrays.toString(freqsValue)				);
//				return posChar;
//			}
//			
//		}
//		System.err.println(d +"\t"+ sumFreq);
//		posChar[1] = GAP;
		return posChar;
//
	}

	
}
