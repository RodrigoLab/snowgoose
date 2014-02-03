package srp.dr.evolution.datatype;

import dr.evolution.datatype.Nucleotides;

public class ShortReads extends Nucleotides {

	
	/**
	 * 
	 */
	private static final long serialVersionUID = 2787317916301857380L;
	public static final String DESCRIPTION = "shortreads";
	public static final int TYPE = 9;
	public static final ShortReads INSTANCE = new ShortReads();

	public static final char SHORTREAD_EMPTY_CHARACTER = '.';
	public static final char SHORTREAD_GAP_CHARACTER = '*';
	/** 
	 * A table to translate state numbers (0-17) into character codes
	 */

	public static final char[] NUCLEOTIDE_CHARS = 
		{ 'A','C','G','T','U','K','M','R','S','W','Y','B','D','H','V','N', 
		UNKNOWN_CHARACTER,GAP_CHARACTER,
		SHORTREAD_EMPTY_CHARACTER, SHORTREAD_GAP_CHARACTER};

	protected ShortReads() {
		stateCount = 4;
		ambiguousStateCount = 18;
	}
	
	@Override
    public char[] getValidChars() {
        return NUCLEOTIDE_CHARS;
    }
	@Override
	public String getDescription() {
		return DESCRIPTION;
	}
	@Override
	public int getType() {
		return TYPE;
	}
	@Override
	public char getChar(int state) {
		return NUCLEOTIDE_CHARS[state];
	}
}
