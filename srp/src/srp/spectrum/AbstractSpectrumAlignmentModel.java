package srp.spectrum;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import srp.dr.evolution.datatype.ShortReads;
import dr.evolution.alignment.SiteList;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.util.Attributable;

//public abstract class AbstractSpectrumAlignmentModel extends AbstractModel implements TaxonList, SiteList{
public abstract class AbstractSpectrumAlignmentModel extends AbstractModel implements TaxonList{
//maybe only can implements TaxonList
	

	private static final long serialVersionUID = -214692337132875593L;
	private static final DataType DATA_TYPE = ShortReads.INSTANCE;
//	private static final int STATE_COUNT = DATA_TYPE.getStateCount();

	public static final char GAP = '-';
	public static final String TAXON_PREFIX = "taxa_";
	
	protected int spectrumLength;

	protected boolean isEdit;//TODO utilise this!!

	public AbstractSpectrumAlignmentModel(String name) {
		super(name);
	}

	public int getSpectrumLength() {
		
		return spectrumLength;
	}
	
	public abstract int getSpectrumCount();
	public abstract AbstractSpectrum getSpectrum(int i);
	public abstract void removeSpectrum(int i);
	public abstract String getSpectrumString(int i);

	
	public double[] getSpecturmFrequencies(int spectrumIndex, int i) {
			return getSpectrum(spectrumIndex).getFrequenciesAt(i);
	}

	public void resetSpectrumOperation() {
		spectrumOperationRecord.setOperation(SpectrumOperation.FULL);
	}

	public SpectrumOperationRecord getSpectrumOperationRecord() {
		return spectrumOperationRecord;
	}

	public SpectrumOperation getSpectrumOperation() {
		return spectrumOperationRecord.getOperation();
	}

	public void setSpectrumOperationRecord(SpectrumOperation op,
			int[] twoSpectrumIndex, int[] swapPositionIndex) {

		spectrumOperationRecord.setRecord(op, twoSpectrumIndex,
				swapPositionIndex);
	}

	public void setSpectrumOperationRecord(SpectrumOperation op, int siteIndex,
			double... delta) {
		spectrumOperationRecord.setRecord(op, siteIndex, delta);
	}

	public void setSpectrumOperationRecord(SpectrumOperation op,
			int spectrumIndex, int siteIndex) {
		spectrumOperationRecord.setRecord(op, spectrumIndex, siteIndex);
	}

	public void setSpectrumOperationRecord(SpectrumOperation op,
			int spectrumIndex, int siteIndex, double delta) {
		spectrumOperationRecord.setRecord(op, spectrumIndex, siteIndex);
	}

	public void setSpectrumOperationRecord(SpectrumOperation op,
			int spectrumIndex, int[] siteIndexs, double... delta) {
		spectrumOperationRecord.setRecord(op, spectrumIndex, siteIndexs, delta);

	}

	public void setSpectrumOperationRecord(SpectrumOperation op,
			int[] spectrumIndexs, int siteIndex) {
		// subcolumn
		spectrumOperationRecord.setRecord(op, spectrumIndexs, siteIndex);

	}

	public void startSpectrumOperation() {
		// System.err.println("\n!!!startSpectrumOperation");
		isEdit = true;
	}

	public void endSpectrumOperation() {
		isEdit = false;
		// System.err.println("!!!!!!!endSpectrumOperation, fireModelCHanged");
		fireModelChanged();
	}
	
	

	@Override
	public void fireModelChanged(){
	//		for (TreeChangedEvent treeChangedEvent : treeChangedEvents) {
			listenerHelper.fireModelChanged(this);//, treeChangedEvent);
	//		}
	//		treeChangedEvents.clear();
		}
	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		System.err.println("Call handleModelChangedEvent in CategorySpectrumAlignmentModel");
		
	}
	@SuppressWarnings("rawtypes")
	@Override
	protected void handleVariableChangedEvent(Variable variable, int index, ChangeType type) {
		System.err.println("Call handleVariableChangedEvent in CategorySpectrumAlignmentModel");
	}
	 
///////////////////////// implementation for  implements TaxonList, SiteList
    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

	@Override
	public int getTaxonCount() {
		return getSpectrumCount();
	}
	

	@Override
	public Taxon getTaxon(int taxonIndex) {
		return getSpectrum(taxonIndex).getTaxon();
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
            throw new IllegalArgumentException("Illegal taxon index:" + taxonIndex);
//            return getSequenceAttribute(taxonIndex, name);
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
	protected SpectrumOperationRecord spectrumOperationRecord;
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

//    @Override
	public DataType getDataType() {
		return DATA_TYPE;
	}
//
//    // **************************************************************
//    // SiteList IMPLEMENTATION
//    // **************************************************************
//
//
//    /**
//     * @return number of sites - getHaplotypeLength()
//     */
//    @Override
//	public int getSiteCount() {
//        return getSpectrumLength();
//    }
//
//    /**
//	 * Gets the pattern index at a particular site
//	 * @param siteIndex
//	 * @return siteIndex, identical to @param
//	 */
//    @Override
//	public int getPatternIndex(int siteIndex) {
//	    return siteIndex;
//	}
//    
//
//    // **************************************************************
//    // PatternList IMPLEMENTATION
//    // **************************************************************
//
//    /**
//     * @return number of patterns
//     */
//    @Override
//	public int getPatternCount() {
//        return getSiteCount();
//    }
// 
//
//    /**
//     * @return number of states for this siteList
//     */
//    @Override
//	public int getStateCount() {
//        return getDataType().getStateCount();
//    }
//
//    /**
//     * Gets the length of the pattern strings which will usually be the
//     * same as the number of taxa
//     *
//     * @return the length of patterns
//     */
//    @Override
//	public int getPatternLength() {
//        return getSpectrumCount();
//    }
//
//    /**
//     * Gets the weight of a site pattern (always 1.0)
//     */
//    @Override
//	public double getPatternWeight(int patternIndex) {
//        return 1.0;
//    }
//
//    /**
//     * @return the array of pattern weights
//     */
//    @Override
//	public double[] getPatternWeights() {
//        double[] weights = new double[getSiteCount()];
//        for (int i = 0; i < weights.length; i++)
//            weights[i] = 1.0;
//        return weights;
//    }
//
//	
//  
//	
    // **************************************************************
    // Alignment IMPLEMENTATION
    // **************************************************************
///	@Override
//	public void setDataType(DataType dataType) {
//		this.dataType = dataType;
//	}
//    
///////////////////////////////////////////    
//    //Methods can't be implemented for spectrum
//	/**
//	 * Gets the pattern as an array of state numbers (one per sequence)
//	 * 
//	 * @return the pattern at patternIndex
//	 */
//	@Override
//	public int[] getPattern(int patternIndex) {
//		throw new IllegalArgumentException(
//				"Not implemented for AbstractSpectrumModel");
//	}
//
//	/**
//	 * @return state at (taxonIndex, patternIndex)
//	 */
//	@Override
//	public int getPatternState(int taxonIndex, int patternIndex) {
//		throw new IllegalArgumentException(
//				"Not implemented for AbstractSpectrumModel");
//	}
//
//	/**
//	 * Gets the pattern of site as an array of state numbers (one per sequence)
//	 * 
//	 * @return the site pattern at siteIndex
//	 */
//	@Override
//	public int[] getSitePattern(int siteIndex) {
//		throw new IllegalArgumentException(
//				"getSitePattern() is not implemented for AbstractSpectrumModel");
//	}
//
//	/**
//	 * @return the frequency of each state
//	 */
//	@Override
//	public double[] getStateFrequencies() {
//		throw new IllegalArgumentException(
//				"Not implemented for AbstractSpectrumModel");
//		//
//		// return PatternList.Utils.empiricalStateFrequencies(this);
//	}
//
//	/**
//	 * @return the sequence state at (taxon, site)
//	 */
//	@Override
//	public int getState(int taxonIndex, int siteIndex) {
//		throw new IllegalArgumentException(
//				"getState() is not implemented for AbstractSpectrumModel");
//	}
///////////////////////// End implementation for  implements TaxonList, SiteList


}
