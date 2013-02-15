package alignment;

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
	
//	HashMap<Integer, String> fullSrp; // seq_id to sequence_with_gaps, aligned to ref, length == ref.length()
//	HashMap<Integer, String> fragmentSrp; // seq_id to sequence_only

	
	int length;
	
	private Integer srpCount;
//	private Pattern p = Pattern.compile("(\\.*)([ACGT*]+)	(\\.*)");
	private Pattern p = Pattern.compile("(\\.*)([ACGT*XYZ]+)(\\.*)");
	
	
	private void init(int length){
		this.length = length;
		
		mapToSrp = new ArrayList[this.length];
		for (int i = 0; i < mapToSrp.length; i++) {
			mapToSrp[i] = new ArrayList<>(); 
		}

		seqNameToSeqID = new HashMap<>();
		shortReads = new ArrayList<>();
//		fullSrp = new HashMap<>();
		srpCount = 0;
		
		
	}
	

	public AlignmentMapping(Alignment srpAlignment){
		init( srpAlignment.getSiteCount() );
		
		for (int i = 0; i < srpAlignment.getSequenceCount(); i++) {
			Sequence s = srpAlignment.getSequence(i);
			addSequence(s);
//					
			
		}
		
		
	}


	public void addSequence(Sequence s) {

		ShortRead srp = new ShortRead(s);
		if (srp.isValid()){
			seqNameToSeqID.put(srp.getName(), srpCount);
			shortReads.add(srpCount, srp);
			
			for (int j = srp.getStart(); j < srp.getEnd(); j++) {
				mapToSrp[j].add(srpCount);
			}
			srpCount++;
		}
		
	}

//
//	private void addSeq(String seqName, String seq) {
//
//		
////		System.out.println(seqName +"\t"+ srpCount);
//		
//		Matcher m = p.matcher(seq);
//
//		if (m.matches() == true){
//		
//			seqNameToSeqID.put(seqName, srpCount);
//			int start = m.start(2);
//			int length = m.end(2) - start;
//			for (int i = m.start(2); i < m.end(2); i++) {
//				mapToSrp[i].add(srpCount);
//			}
//			
//			fullSrp.put(srpCount, seq);
//			srpCount++;
////			System.out.println(start +"\t"+  end + seq.charAt(start) +"\t"+ seq.charAt(end));			
////			System.out.println(m.start(2) +"\t"+  m.end(2));
////			String srp2 = seq.substring(m.start(2), m.end(2));
////			System.out.println(srp.equalsIgnoreCase(srp2));
////			System.out.println(m.group(3));
////			System.out.println(m.toMatchResult());
//
//		}
//
//	}
	public ShortRead getShortRead(int i){
		return shortReads.get(i);
	}
	
	
	public char getShortReadCharAt(int index, int c){
		return shortReads.get(index).getFullSrpCharAt(c);
	}
	
	public String toString() {
		String s = "";
		for (int i = 0; i < mapToSrp.length; i++) {
			s +=mapToSrp[i].toString();
			s += "\n";
			
		}
		return s;
		
	}


	public int getSrpCount() {
		return srpCount;
	}


	public String getFragmentSrp(int i) {
		
		return shortReads.get(i).getFragmentSrp();
	}

	public String getFullSrp(int i) {
		
		return shortReads.get(i).getFullSrp();
	}


	public int getStart(int i) {
	
		return shortReads.get(i).getStart();
	}


	public int getEnd(int i) {

		return shortReads.get(i).getEnd();
	}


	public int getSrpLength(int i) {
		return shortReads.get(i).getLength();
	}
	public int getLength(){
		return length;
	}
}
