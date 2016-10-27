package srp.evolution;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import srp.dr.evolution.datatype.ShortReads;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.util.Attributable;

public abstract class AbstractAlignmentModel extends AbstractModel implements TaxonList, Attributable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 2358008010540023190L;
	
	public static final DataType DATA_TYPE = ShortReads.INSTANCE;
	public static final char GAP = DataType.GAP_CHARACTER;
	public static final String TAXON_PREFIX = "hap_";
	
	protected boolean isEdit;
	protected DataType dataType;
	protected OperationRecord operationRecord;
	
	protected String id = null;
	protected Attributable.AttributeHelper attributes = null;
	
	protected int sequenceCount;
	protected int sequenceLength;

	public AbstractAlignmentModel(String name) {
		super(name);
		operationRecord = new OperationRecord();
		dataType = DATA_TYPE;
	}

	public int getSequenceCount() {
		return sequenceCount;
	}
	public int getSequenceLength(){
		return sequenceLength;
	}
	
		
	
	
	public void startAlignmentModelOperation() {
		// System.err.println("\n!!!startSpectrumOperation");
		isEdit = true;
	}

	public void endAlignmentModelOperation() {
		isEdit = false;
		// System.err.println("!!!!!!!endSpectrumOperation, fireModelCHanged");
		fireModelChanged();
	}

	public void setDataType(DataType dataType) {
		this.dataType = dataType;
	}
	
	public DataType getDataType() {
		return dataType;
	}


	public void resetOperation() {
		operationRecord.setOperation(OperationType.FULL);
	}

	public OperationRecord getOperationRecord() {
		return operationRecord;
	}

	public OperationType getOperation() {
		return operationRecord.getOperation();
	}

	public void setOperationRecord(OperationType op,
			int[] twoSpectrumIndex, int[] swapPositionIndex) {
		//Recom
		operationRecord.setRecord(op, twoSpectrumIndex, swapPositionIndex);
	}

	public void setOperationRecord(OperationType op,
			int spectrumIndex, int siteIndex) {
		//Single
		operationRecord.setRecord(op, spectrumIndex, siteIndex);
	}

	public void setOperationRecord(OperationType op,
			int spectrumIndex, int[] siteIndexs) {
		//Multi
		operationRecord.setRecord(op, spectrumIndex, siteIndexs);

	}

	public void setOperationRecord(OperationType op,
			int[] spectrumIndexs, int siteIndex) {
		//Column, subcolumn
		operationRecord.setRecord(op, spectrumIndexs, siteIndex);

	}
	@Deprecated
	public void setOperationRecord(OperationType op, int siteIndex,
			double... delta) {
		operationRecord.setRecord(op, siteIndex, delta);
	}

	@Deprecated
	public void setOperationRecord(OperationType op,
			int spectrumIndex, int siteIndex, double delta) {
		operationRecord.setRecord(op, spectrumIndex, siteIndex);
	}
	@Deprecated
	public void setOperationRecord(OperationType op,
			int spectrumIndex, int[] siteIndexs, double... delta) {
		operationRecord.setRecord(op, spectrumIndex, siteIndexs, delta);

	}

	@Override
	public void fireModelChanged() {
		// for (TreeChangedEvent treeChangedEvent : treeChangedEvents) {
		listenerHelper.fireModelChanged(this);// , treeChangedEvent);
		// }
		// treeChangedEvents.clear();
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		System.err
				.println("Call handleModelChangedEvent in CategorySpectrumAlignmentModel");

	}

	@SuppressWarnings("rawtypes")
	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		System.err
				.println("Call handleVariableChangedEvent in CategorySpectrumAlignmentModel");
	}

	@Override
	protected void acceptState() {
		//Do nothing
	}
    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************


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



	/**
	 * Sets an named attribute for this object.
	 *
	 * @param name  the name of the attribute.
	 * @param value the new value of the attribute.
	 */
	@Override
	public void setAttribute(String name, Object value) {
	    if (attributes == null)
	        attributes = new Attributable.AttributeHelper();
	    attributes.setAttribute(name, value);
	}

	/**
	 * @param name the name of the attribute of interest.
	 * @return an object representing the named attributed for this object.
	 */
	@Override
	public Object getAttribute(String name) {
	    if (attributes == null)
	        return null;
	    else
	        return attributes.getAttribute(name);
	}

	/**
	 * @return an iterator of the attributes that this object has.
	 */
	@Override
	public Iterator<String> getAttributeNames() {
	    if (attributes == null)
	        return null;
	    else
	        return attributes.getAttributeNames();
	}


}