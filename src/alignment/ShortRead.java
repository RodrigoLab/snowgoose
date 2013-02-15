package alignment;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import dr.evolution.sequence.Sequence;

public class ShortRead {

	private String fullSrp;
	private String fragmentSrp;
	private String srpName;
	
	private int start;
	private int end;
	
	private int length;
	
	private int srpID;
	
	private Pattern p = Pattern.compile("(\\.*)([ACGT*XYZ]+)(\\.*)");
//	private Pattern p = Pattern.compile("(\\.*)([ACGT*]+)(\\.*)");
	//TODO: change back later
	private boolean isValid = false;
	
	
	public ShortRead(Sequence seq) {
		this.fullSrp = seq.getSequenceString();
		this.srpName = seq.getTaxon().getId();

		Matcher m = p.matcher(fullSrp);

		if (m.matches() == true){
			isValid = true;
			setStartEnd(m.start(2), m.end(2));
			fragmentSrp = fullSrp.substring(start,end);

//			System.out.println(fragmentSrp.equals(m.group(2)));
//			System.out.println(fragmentSrp.length() +"\t"+ (end-start));
					

		}

		
		
	}
	
	public ShortRead(String fullSrp, String fragmentSrp, String name,
			int srpID, int start, int length, int end) {
		super();
		this.fullSrp = fullSrp;
		this.fragmentSrp = fragmentSrp;
		this.srpName = name;
		this.srpID = srpID;
		this.start = start;
		this.length = length;
		this.end = end;
	}

	private void setStartEnd(int start, int end) {
		this.start = start;
		this.end = end;
		this.length = this.end - this.start;
		
	}
	
	public char getFullSrpCharAt(int i){
		return fullSrp.charAt(i);
	}

	public boolean isValid() {
		return isValid ;
	}
	public String getName(){
		return srpName;
	}
	public int getStart(){
		return start;
	}
	public int getEnd(){
		return end;
	}

	public String getFullSrp() {
		return fullSrp;
	}
	public String getFragmentSrp(){
		return fragmentSrp;
	}

	public int getLength() {
		return length;
	}
	
}