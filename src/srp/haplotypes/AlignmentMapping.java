package srp.haplotypes;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.math.MathUtils;

public class AlignmentMapping {

	private static final int GAP = '-';
	private static final char[] VALID_CHARS = Nucleotides.INSTANCE.getValidChars();

	private ArrayList<Integer>[] mapToSrp; // each [] = position, each ArrayList
											// = map to which read
	private HashMap<String, Integer> seqNameToSeqID; // map sequence_name >xxxto int

	private ArrayList<Character>[] listOfAvailableChar;
	private ArrayList<ShortRead> shortReads;

	private int haplotypeLength;
	private Integer srpCount;
	private double[] cumFreq;
	
	
	private void init(int l) {
		haplotypeLength = l;
		srpCount = 0;

		cumFreq = new double[haplotypeLength];

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
		for (int i = 0; i < setsOfAvailableChar.length; i++) {
			setsOfAvailableChar[i] = new HashSet<Character>();
		}

		for (int i = 0; i < srpAlignment.getSequenceCount(); i++) {
			
			Sequence s = srpAlignment.getSequence(i);
			addSequence(s, setsOfAvailableChar);
		}
		
		for (int i = 0; i < haplotypeLength; i++) {
			listOfAvailableChar[i] = new ArrayList<Character>(setsOfAvailableChar[i]);
		}
		
		cumFreq = new double[] { cumFreq['A'], cumFreq['C'],
				cumFreq['G'], cumFreq['T'] };


		for (int i = 1; i < cumFreq.length; i++) {
			cumFreq[i] = cumFreq[i] + cumFreq[i-1];  
		}
		double sum = cumFreq[3];
		for (int i = 0; i < cumFreq.length; i++) {
			cumFreq[i] /= sum;  
		}
	}


	private void addSequence(Sequence s,
			HashSet<Character>[] setsOfAvailableChar) {

		ShortRead srp = new ShortRead(s);
		if (srp.getIsValid()){
			seqNameToSeqID.put(srp.getName(), srpCount);
			shortReads.add(srpCount, srp);
			
			for (int j = srp.getStart(); j < srp.getEnd(); j++) {
				mapToSrp[j].add(srpCount);
				char c = srp.getFullSrpCharAt(j);
				if (c!= GAP){
					setsOfAvailableChar[j].add(c);
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
	
	public int[] getNextBaseUniform() {
		int newChar = GAP ;
		int pos = MathUtils.nextInt(haplotypeLength);
		int size = listOfAvailableChar[pos].size();
		if (size != 0) {
			newChar = listOfAvailableChar[pos].get(MathUtils.nextInt(size));
		}
		
		return new int[]{pos, newChar};
	}


	public int[] getNextBaseEmpirical() {

		int newChar = GAP;
		int pos = MathUtils.nextInt(haplotypeLength);
		
		double d = MathUtils.nextDouble();
		for (int i = 0; i < cumFreq.length; i++) {
			if (d <= cumFreq[i]) {
				newChar = VALID_CHARS[i];
				break;
			}
			
		}
		

		return new int[] { pos, newChar };
	}
}
