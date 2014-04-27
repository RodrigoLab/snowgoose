package srp.evolution.haplotypes;

import java.util.ArrayList;

import srp.evolution.AbstractAlignmentModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

public abstract class AbstractHaplotypeModelArrayList extends AbstractAlignmentModel implements Alignment {

	
	public AbstractHaplotypeModelArrayList(String name) {
		super(name);
	}

	private static final long serialVersionUID = -214692337132875593L;


	protected int haplotypeLength;
	protected ArrayList<Haplotype> haplotypes;

	public String getHaplotypeString(int i) {
		return haplotypes.get(i).getSequenceString();
		
	}

	public int getHaplotypeLength() {
	
		return haplotypeLength;
	}

	public int getHaplotypeCount() {
		return haplotypes.size();
	}

	public Haplotype getHaplotype(int i){
		return haplotypes.get(i);
	}
	
	public void removeHaplotype(int i) {
		haplotypes.remove(i);
		
	}


	public char getHaplotypeCharAt(int hapIndex, int charIndex) {
		return haplotypes.get(hapIndex).getChar(charIndex);
	}
	
    // **************************************************************
    // SequenceList IMPLEMENTATION
    // **************************************************************
	
	/*
	 * Call getHaplotypesCount(), return haplotypesCount
	 */
	@Override
	public int getSequenceCount() {
		return getHaplotypeCount();
	}

	@Override
	public Sequence getSequence(int i) {
		return getHaplotype(i);
	}

	
	@Override
	public void setSequenceAttribute(int index, String name, Object value) {
		Sequence sequence = getSequence(index);
        sequence.setAttribute(name, value);
	}

	
	@Override
	public Object getSequenceAttribute(int index, String name) {
		Sequence sequence = getSequence(index);
        return sequence.getAttribute(name);
	}

    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

	@Override
	public int getTaxonCount() {
		return getHaplotypeCount();
	}
	
	@Override
	public Taxon getTaxon(int taxonIndex) {
		return getHaplotype(taxonIndex).getTaxon();
	}



    // **************************************************************
    // SiteList IMPLEMENTATION
    // **************************************************************


    /**
     * @return number of sites - getHaplotypeLength()
     */
    @Override
	public int getSiteCount() {
        return getHaplotypeLength();
    }

    /**
     * Gets the pattern of site as an array of state numbers (one per sequence)
     *
     * @return the site pattern at siteIndex
     */
    @Override
	public int[] getSitePattern(int siteIndex) {

    	int n = getHaplotypeCount();

    	int[] pattern = new int[n];
        for (int i = 0; i < n; i++) {
            Haplotype hap = getHaplotype(i);

            if (siteIndex >= hap.getLength())
                pattern[i] = dataType.getGapState();
            else
                pattern[i] = hap.getState(siteIndex);
        }

        return pattern;
    }
	/**
	 * Gets the pattern index at a particular site
	 * @param siteIndex
	 * @return siteIndex, identical to @param
	 */
    @Override
	public int getPatternIndex(int siteIndex) {
	    return siteIndex;
	}
    
    /**
     * @return the sequence state at (taxon, site)
     */
    @Override
	public int getState(int taxonIndex, int siteIndex) {
        Haplotype hap = getHaplotype(taxonIndex);

        if (siteIndex >= hap.getLength()) {
            return dataType.getGapState();
        }

        return hap.getState(siteIndex);
    }

//    /**
//     */
//    public void setState(int taxonIndex, int siteIndex, int state) {
//
//        Sequence seq = getSequence(taxonIndex);
//
//        if (siteIndex >= seq.getLength()) {
//            throw new IllegalArgumentException();
//        }
//
//        seq.setState(siteIndex, state);
//    }

    // **************************************************************
    // PatternList IMPLEMENTATION
    // **************************************************************

    /**
     * @return number of patterns
     */
    @Override
	public int getPatternCount() {
        return getSiteCount();
    }
 

    /**
     * @return number of states for this siteList
     */
    @Override
	public int getStateCount() {
        return getDataType().getStateCount();
    }

    /**
     * Gets the length of the pattern strings which will usually be the
     * same as the number of taxa
     *
     * @return the length of patterns
     */
    @Override
	public int getPatternLength() {
        return getSequenceCount();
    }

    /**
     * Gets the pattern as an array of state numbers (one per sequence)
     *
     * @return the pattern at patternIndex
     */
    @Override
	public int[] getPattern(int patternIndex) {

        return getSitePattern(patternIndex);
    }

    /**
     * @return state at (taxonIndex, patternIndex)
     */
    @Override
	public int getPatternState(int taxonIndex, int patternIndex) {
    	return getState(taxonIndex, patternIndex);
    }

    /**
     * Gets the weight of a site pattern (always 1.0)
     */
    @Override
	public double getPatternWeight(int patternIndex) {
        return 1.0;
    }

    /**
     * @return the array of pattern weights
     */
    @Override
	public double[] getPatternWeights() {
        double[] weights = new double[getSiteCount()];
        for (int i = 0; i < weights.length; i++)
            weights[i] = 1.0;
        return weights;
    }

	
    /**
     * @return the frequency of each state
     */
    @Override
	public double[] getStateFrequencies() {
        return PatternList.Utils.empiricalStateFrequencies(this);
    }


	
    // **************************************************************
    // Alignment IMPLEMENTATION
    // **************************************************************



	
	/**
	 * call getHaplotypeString(sequenceIndex);
	 */
	@Override
	public String getAlignedSequenceString(int sequenceIndex) {
		
		return getHaplotypeString(sequenceIndex);
	}

	/**
	 * call getHaplotypeString(sequenceIndex);
	 */
	@Override
	public String getUnalignedSequenceString(int sequenceIndex) {
		
		return getHaplotypeString(sequenceIndex);
	}


}
