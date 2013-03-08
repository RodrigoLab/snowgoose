package srp.haplotypes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import dr.evolution.alignment.Alignment;
import dr.evolution.sequence.Sequence;

public class AlignmentMapping {

	ArrayList<Integer>[] mapToSrp; // each [] = position, each ArrayList = map to which read
	HashMap<String, Integer> seqNameToSeqID; // map sequence_name >xxx to int

	ArrayList<ShortRead> shortReads;

	private int length;
	private Integer srpCount;

	
	private void init(int length){
		this.length = length;
		
		mapToSrp = new ArrayList[this.length];
		for (int i = 0; i < mapToSrp.length; i++) {
			mapToSrp[i] = new ArrayList<Integer>(); 
		}

		seqNameToSeqID = new HashMap<String, Integer>();
		shortReads = new ArrayList<ShortRead>();
//		fullSrp = new HashMap<>();
		srpCount = 0;
		
		
	}
	

	public AlignmentMapping(Alignment srpAlignment){
		init( srpAlignment.getSiteCount() );
		
		for (int i = 0; i < srpAlignment.getSequenceCount(); i++) {
			Sequence s = srpAlignment.getSequence(i);
			addSequence(s);
		}
	}


	public void addSequence(Sequence s) {

		ShortRead srp = new ShortRead(s);
		if (srp.getIsValid()){
			seqNameToSeqID.put(srp.getName(), srpCount);
			shortReads.add(srpCount, srp);
			
			for (int j = srp.getStart(); j < srp.getEnd(); j++) {
				mapToSrp[j].add(srpCount);
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
		return length;
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

}
