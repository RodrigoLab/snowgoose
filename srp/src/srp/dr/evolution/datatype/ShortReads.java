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
	 * A table to translate state numbers (0-19) into character codes
	 */

			
	public static final char[] NUCLEOTIDE_CHARS = 
		{ 'A','C','G','T','U','K','M','R','S','W','Y','B','D','H','V','N', 
		UNKNOWN_CHARACTER,GAP_CHARACTER,
		SHORTREAD_EMPTY_CHARACTER, SHORTREAD_GAP_CHARACTER};

	
	/** 
	 * This table maps nucleotide characters into state codes (0-19)
	 * Nucleotides go ACGTURYMWSKBDHVN?-", Other letters are mapped to ?.
	 * ? and - are mapped to themselves. All other chars are mapped to -.
	 *
	 	?	63
		-	45
		.	46	SHORTREAD_EMPTY_CHARACTER
		*	42	SHORTREAD_GAP_CHARACTER
	 */
	public static final int NUCLEOTIDE_STATES[] = {
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 0-15
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 16-31
	//                                *         -  .
		17,17,17,17,17,17,17,17,17,17,19,17,17,17,18,17,	// 32-47
	//                                                ?
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,16,	// 48-63
	//	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
		17, 0,11, 1,12,16,16, 2,13,16,16,10,16, 7,15,16,	// 64-79
	//	 p  q  R  S  T  U  V  W  x  Y  z
		16,16, 5, 9, 3, 3,14, 8,16, 6,16,17,17,17,17,17,	// 80-95
	//	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
		17, 0,11, 1,12,16,16, 2,13,16,16,10,16, 7,15,16,	// 96-111
	//	 p  q  R  S  T  U  V  W  x  Y  z
		16,16, 5, 9, 3, 3,14, 8,16, 6,16,17,17,17,17,17		// 112-127
	};

	
	
	protected ShortReads() {
		stateCount = 4;
		ambiguousStateCount = 19;
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
	@Override
	public int getState(char c) {
		return NUCLEOTIDE_STATES[c];
	}
}
