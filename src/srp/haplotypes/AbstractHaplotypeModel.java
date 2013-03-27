package srp.haplotypes;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import dr.inference.model.AbstractModel;
import dr.util.Attributable;

public abstract class AbstractHaplotypeModel  extends AbstractModel implements Alignment {

	
	public AbstractHaplotypeModel(String name) {
		super(name);
	}

	private static final long serialVersionUID = -214692337132875593L;
	private DataType dataType = null;

	protected int haplotypesLength;
	protected ArrayList<Haplotype> haplotypes;

	public String getHaplotypeString(int i) {
	
		return getHaplotype(i).getSequenceString();
		
	}

	public int getHaplotypeLength() {
	
		return haplotypesLength;
	}

	public int getHaplotypeCount() {
		return haplotypes.size();
	}

	public Haplotype getHaplotype(int i){
		return haplotypes.get(i);
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

	/**
     * @return the ID of the taxon of the ith sequence. If it doesn't have
     *         a taxon, returns the ID of the sequence itself.
     */
    @Override
	public String getTaxonId(int taxonIndex) {
        Taxon taxon = getTaxon(taxonIndex);
        if (taxon != null)
            return taxon.getId();
        else
            throw new IllegalArgumentException("Illegal taxon index:" + taxonIndex);
    }

    /**
     * returns the index of the taxon with the given id.
     */
    @Override
	public int getTaxonIndex(String id) {
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            if (getTaxonId(i).equals(id)) return i;
        }
        return -1;
    }

    /**
     * returns the index of the given taxon.
     * must be the same object
     */
    @Override
	public int getTaxonIndex(Taxon taxon) {
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            if (getTaxon(i) == taxon) return i;
        }
        return -1;
    }

    @Override
	public List<Taxon> asList() {
        List<Taxon> taxa = new ArrayList<Taxon>();
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            taxa.add(getTaxon(i));
        }
        return taxa;
    }

	
//    /**
//     * Sets an named attribute for the taxon of a given sequence. If the sequence
//     * doesn't have a taxon then the attribute is added to the sequence itself.
//     *
//     * @param taxonIndex the index of the taxon whose attribute is being set.
//     * @param name       the name of the attribute.
//     * @param value      the new value of the attribute.
//     */
//    public void setTaxonAttribute(int taxonIndex, String name, Object value) {
//        Taxon taxon = getTaxon(taxonIndex);
//        if (taxon != null)
//            taxon.setAttribute(name, value);
//        else
//            setSequenceAttribute(taxonIndex, name, value);
//    }

    /**
     * @param taxonIndex the index of the taxon whose attribute is being fetched.
     * @param name       the name of the attribute of interest.
     * @return an object representing the named attributed for the given taxon.
     */
    @Override
	public Object getTaxonAttribute(int taxonIndex, String name) {
    	Taxon taxon = getTaxon(taxonIndex);
        if (taxon != null)
            return taxon.getAttribute(name);
        else
            return getSequenceAttribute(taxonIndex, name);
    }

    @Override
	public Iterator<Taxon> iterator() {
        return new Iterator<Taxon>() {
            private int index = -1;

            @Override
			public boolean hasNext() {
                return index < getTaxonCount() - 1;
            }

            @Override
			public Taxon next() {
                index++;
                return getTaxon(index);
            }

            @Override
			public void remove() { /* do nothing */ }
        };
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

	@Override
	public DataType getDataType() {
		return dataType;
	}
	
    // **************************************************************
    // Alignment IMPLEMENTATION
    // **************************************************************

	@Override
	public void setDataType(DataType dataType) {
		this.dataType = dataType;
	}

	
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
	
    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

	protected String id = null;

	/**
	 * @return the id.
	 */
	@Override
	public String getId() {
		return id;
	}

	/**
	 * Sets the id.
	 */
	@Override
	public void setId(String id) {
		this.id = id;
	}


    // **************************************************************
    // Attributable IMPLEMENTATION
    // **************************************************************
	private Attributable.AttributeHelper attributes = null;
    /**
     * Sets an named attribute for this object.
     *
     * @param name  the name of the attribute.
     * @param value the new value of the attribute.
     */
    public void setAttribute(String name, Object value) {
        if (attributes == null)
            attributes = new Attributable.AttributeHelper();
        attributes.setAttribute(name, value);
    }

    /**
     * @param name the name of the attribute of interest.
     * @return an object representing the named attributed for this object.
     */
    public Object getAttribute(String name) {
        if (attributes == null)
            return null;
        else
            return attributes.getAttribute(name);
    }

    /**
     * @return an iterator of the attributes that this object has.
     */
    public Iterator<String> getAttributeNames() {
        if (attributes == null)
            return null;
        else
            return attributes.getAttributeNames();
    }

}
