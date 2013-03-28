package srp.haplotypes;

import dr.evolution.datatype.DataType;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

public class Haplotype extends Sequence {

	/**
	 * 
	 */
	private static final long serialVersionUID = -1481408524381652240L;
	StringBuilder sequenceString;

	public Haplotype() {
		sequenceString = new StringBuilder();
	}

	public Haplotype(String sequence) {
		this();
		setSequenceString(sequence);
	}

	public Haplotype(Sequence sequence) {
		this(sequence.getTaxon(), sequence.getSequenceString());
	}

	public Haplotype(Taxon taxon, String sequence) {
		this();
		setTaxon(taxon);
		setSequenceString(sequence);
	}

	
	public void setCharAt(int index, char newChar) {
		sequenceString.setCharAt(index, newChar);
	}

	public void setHaplotypeString(String haplotypeString) {
		setSequenceString(haplotypeString);

	}
	
	public char charAt(int index){
		return sequenceString.charAt(index);
	}
	
	// **************************************
	// OVERRIDE ALL (almost all) methods
	// Do NOT call setState()!!
	// ************************************
	

    /**
     * @return the DataType of the sequences.
     */
    @Override
	public DataType getDataType() {
        return dataType;
    }

    /**
     * @return the length of the sequences.
     */
    @Override
	public int getLength() {
        return sequenceString.length();
    }

    /**
     * @return a String containing the sequences.
     */
    @Override
	public String getSequenceString() {
        return sequenceString.toString();
    }

    /**
     * @return a char containing the state at index.
     */
    @Override
	public char getChar(int index) {
        return charAt(index);
    }

    /**
     * @return the state at site index.
     */
    @Override
	public int getState(int index) {
        return dataType.getState(sequenceString.charAt(index));
    }

    /**
     */
//	public void setState(int index, int state) {
//
//        sequenceString.setCharAt(index, dataType.getChar(state));
//    }

    /**
     * Characters are copied from the sequences into the destination character array dst.
     */
    @Override
	public void getChars(int srcBegin, int srcEnd, char[] dst, int dstBegin) {
        sequenceString.getChars(srcBegin, srcEnd, dst, dstBegin);
    }

    /**
     * Set the DataType of the sequences.
     */
    @Override
	public DataType guessDataType() {
        return DataType.guessDataType(sequenceString.toString());
    }

    /**
     * Set the sequences using a string.
     */
    @Override
	public void setSequenceString(String sequence) {
        sequenceString.setLength(0);
        sequenceString.append(sequence.toUpperCase());
    }

    /**
     * Append a string to the sequences.
     */
    @Override
	public void appendSequenceString(String sequence) {
        sequenceString.append(sequence);
    }

    /**
     * Insert a string into the sequences.
     */
    @Override
	public void insertSequenceString(int offset, String sequence) {
        sequenceString.insert(offset, sequence);
    }




}