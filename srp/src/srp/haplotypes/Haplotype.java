package srp.haplotypes;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;

public class Haplotype extends Sequence {

	/**
	 * 
	 */
	private static final long serialVersionUID = -1481408524381652240L;
	
//	@Deprecated 
//	StringBuilder sequenceString;

	
	char[] haplotype;
	char[] storedHaplotype;
	int haplotypeLength;

	public static final DataType DATATYPE = Nucleotides.INSTANCE;
	public static final int STATE_COUNT = DATATYPE.getStateCount();

	public Haplotype(int haplotypeLength) {
		setDataType(DATATYPE);
		sequenceString = null;
		this.haplotypeLength = haplotypeLength;
		haplotype = new char[this.haplotypeLength];
		storedHaplotype = new char[this.haplotypeLength];
		initHaplotype();
	}

	private void initHaplotype() {

		for (int i = 0; i < haplotypeLength; i++) {
			int state = MathUtils.nextInt(STATE_COUNT);
			char newChar = DATATYPE.getChar(state);
			haplotype[i] = newChar;
		}
		storeState();
		
	}

	public Haplotype(String sequence) {
		this(sequence.length());
		setSequenceString(sequence.toUpperCase());
		storeState();
	}

	public Haplotype(Sequence sequence) {
		this(sequence.getTaxon(), sequence.getSequenceString());
	}

	public Haplotype(Taxon taxon, String sequence) {
		this(sequence);
		setTaxon(taxon);
	}

	public Haplotype(Taxon taxon, int haplotypeLength) {
		this(haplotypeLength);
		setTaxon(taxon);
	}

	public void setCharAt(int index, int newChar) {
		setCharAt(index, (char) newChar); 
	}
	
	public void setCharAt(int index, char newChar) {
		haplotype[index]=newChar;
	}

	@Deprecated
	public char replaceCharAt(int index, int newChar){
		char oldChar = getChar(index);
		setCharAt(index, (char) newChar);
		return oldChar;
	}
	
//	@Override
	public char getStoredChar(int index) {
        return storedHaplotype[index];
    }
	public int getStoredState(int index) {
		return dataType.getState(storedHaplotype[index]);
//        return storedHaplotype[index];
    }
	// **************************************
	// OVERRIDE ALL (almost all) methods
	// Do NOT call setState()!!
	// ************************************
	
	@Override
    public void setDataType(DataType dataType) {
        this.dataType = DATATYPE;
    }

    /**
     * @return the length of the sequences.
     */
    @Override
	public int getLength() {
        return haplotypeLength;
    }

    /**
     * @return a String containing the sequences.
     */
    @Override
	public String getSequenceString() {
    	return String.valueOf(haplotype);
    	
    }

    /**
     * @return a char containing the state at index.
     */
    @Override
	public char getChar(int index) {
        return haplotype[index];
    }

    /**
     * @return the state at site index.
     */
    @Override
	public int getState(int index) {
        return dataType.getState(haplotype[index]);
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
        System.arraycopy(haplotype, srcBegin, dst, dstBegin, srcEnd - srcBegin);
    }

    /**
     * Set the DataType of the sequences.
     */
    @Override
	public DataType guessDataType() {
        DataType guessDataType = DataType.guessDataType(String.valueOf(haplotype));
        if(guessDataType.getName().equals(DATATYPE.getName())){
        	return DATATYPE;
        }
        else{
        	throw new IllegalArgumentException("Only support "+DATATYPE.getName()+". Please check your haplotypes");
        }
        
    }

    /**
     * Set the sequences using a string.
     */
    @Override
	public void setSequenceString(String sequence) {
    	if(sequence.length() != haplotypeLength){
			throw new IllegalArgumentException("Invalid sequence length: "
					+ sequence.length()
					+ ". Haplotype length must be equal to " + haplotypeLength);
		}
    	else{
    		System.arraycopy(sequence.toCharArray(), 0, haplotype, 0, haplotypeLength);
    	}
    }

    /**
     * Append a string to the sequences.
     */
    @Override
	public void appendSequenceString(String sequence) {
        throw new IllegalArgumentException("Can not alter haplotype length");
    }

    /**
     * Insert a string into the sequences.
     */
    @Override
	public void insertSequenceString(int offset, String sequence) {
    	throw new IllegalArgumentException("Can not alter haplotype length");
    }

	public void storeState(int index) {
		storedHaplotype[index] = haplotype[index];
	}
	public void restoreState(int index) {
		haplotype[index] = storedHaplotype[index];
	}
	
	public void storeState() {
		System.arraycopy(haplotype, 0, storedHaplotype, 0, haplotypeLength);
	}

	public void restoreState() {
		char[] temp = storedHaplotype;
		storedHaplotype = haplotype;
		haplotype = temp;
	}

	public static Haplotype duplicateHaplotype(Haplotype oldHaplotype) {
		
		Taxon newTaxon = oldHaplotype.getTaxon();
		Haplotype newHaplotype = new Haplotype(newTaxon, oldHaplotype.getSequenceString());

		return newHaplotype;
	}

	public String getStoredSequenceString() {
    	return String.valueOf(storedHaplotype);
	}



}