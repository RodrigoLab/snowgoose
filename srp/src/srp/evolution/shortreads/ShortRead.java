package srp.evolution.shortreads;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import dr.evolution.sequence.Sequence;

public class ShortRead {

	private static Pattern p = Pattern.compile("(\\.*)([ACGT*]+)(\\.*)");
	
	private String fullSrp;
	private String fragmentSrp;
	private String srpName;

	private int start;
	private int end;
	private int length;

	private boolean isValid = false;

	private Integer lengthInteger;

//	private char[] fullSrp2;

	public ShortRead(Sequence seq) {
		this.fullSrp = seq.getSequenceString();
		this.srpName = seq.getTaxon().getId();

		Matcher m = p.matcher(fullSrp);

		if (m.matches() == true) {
			isValid = true;
			setStartEnd(m.start(2), m.end(2));
			fullSrp = fullSrp.replaceAll("\\*", "-");
			fragmentSrp = fullSrp.substring(start, end);
			
			
			// System.out.println(fragmentSrp.equals(m.group(2)));
			// System.out.println(fragmentSrp.length() +"\t"+ (end-start));

		}
//		fullSrp2 = fullSrp.toCharArray();
	}

	private void setStartEnd(int start, int end) {
		this.start = start;
		this.end = end;
		this.length = this.end - this.start;
		lengthInteger = length;
	}

	public char getFullSrpCharAt(int i) {
		return fullSrp.charAt(i);
//		return fullSrp2[i];
	}

	public boolean getIsValid() {
		return isValid;
	}

	public String getName() {
		return srpName;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public String getFullSrp() {
		return fullSrp;
	}

	public String getFragmentSrp() {
		return fragmentSrp;
	}

	public int getLength() {
		return length;
	}

	public Integer getLengthInteger() {

		return lengthInteger;
		
	}

}