package srp.spectrum;

import java.util.Iterator;

import srp.dr.evolution.datatype.ShortReads;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.util.Attributable;

public abstract class AbstractSpectrum extends AbstractModel implements Attributable  {
//public abstract class AbstractSpectrum extends AbstractModel  {

	/**
	 * 
	 */
	private static final long serialVersionUID = -443784564247356367L;
	protected static final DataType DATA_TYPE = ShortReads.INSTANCE;
	protected static final int STATE_COUNT = DATA_TYPE.getStateCount();

	protected Taxon taxon;
	protected int storeSiteIndex;
    // **************************************************************
    // Attributable IMPLEMENTATION
    // **************************************************************

    private Attributable.AttributeHelper attributes = null;

	public AbstractSpectrum(String name) {
		super(name);
	}

	public abstract int getLength();

	/**
	 * Sets a taxon for this sequences.
	 *
	 * @param taxon the taxon.
	 */
	public void setTaxon(Taxon taxon) {
	    this.taxon = taxon;
	}

	/**
	 * @return the taxon for this sequences.
	 */
	public Taxon getTaxon() {
	
	    return taxon;
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
	// **************************************************************
	// AbstractModel IMPLEMENTATION
	// **************************************************************

	@Override
	protected void handleModelChangedEvent(Model model, Object object,
			int index) {
				System.err.println("Call handleModelChangedEvent");
			}

	@SuppressWarnings("rawtypes")
	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
			//		System.err.println("Call handleVariableChangedEvent Specturm\t"+variable.getVariableName() +"\t"+ index);
					//TODO implement when secptra changed
				}

	public void setStoreSiteIndex(int s) {
	//		System.err.println("setStoreSiteIndex\t"+s);
			storeSiteIndex = s;
		}

	@Override
	protected void storeState() {
	//		System.err.println("storeState Specturm: "+storeSiteIndex);
	
			getSpectra(storeSiteIndex).storeValues();
	//		for (int i = 0; i < spectrum.size(); i++) {
	//			storeSpectrum.set(i, spectrum.get(i));
	//		}
		}

	@Override
	protected void restoreState() {
	//		System.err.println("restoreState Spectrum: "+storeSiteIndex);
			getSpectra(storeSiteIndex).restoreValues();
	
	//		for (int i = 0; i < spectrum.size(); i++) {
	//			spectrum.set(i, storeSpectrum.get(i));
	//		}
	
		}

	@Override
	protected void acceptState() {
		// Do nothing
		
	}
	
	abstract AbstractSpectra getSpectra(int i);

}