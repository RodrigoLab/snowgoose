package srp.spectrum;

import java.text.NumberFormat;
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

	private static final long serialVersionUID = -443784564247356367L;
	protected static final DataType DATA_TYPE = ShortReads.INSTANCE;
	protected static final int STATE_COUNT = DATA_TYPE.getStateCount();

	protected Taxon taxon;
	protected int storeSiteIndex;
	protected int spectrumLength;
	
	public abstract AbstractSpectra getSpectra(int i);

	public AbstractSpectrum(String name) {
		super(name);
	}
	
	public int getLength(){
		return spectrumLength;
	}
	
	public double getFrequency(int site, int state) {
	    return getSpectra(site).getFrequency(state);
	}

	public double[] getFrequenciesAt(int site) {
		return getSpectra(site).getFrequencies();
	}

	public void setTaxon(Taxon taxon) {
	    this.taxon = taxon;
	}

	/**
	 * @return the taxon for this sequences.
	 */
	public Taxon getTaxon() {
	    return taxon;
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		System.err.println("Call handleModelChangedEvent");
	}

	@SuppressWarnings("rawtypes")
	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		// System.err.println("Call handleVariableChangedEvent Specturm\t"+variable.getVariableName()
		// +"\t"+ index);
		// TODO implement when secptra changed
	}

	public void setStoreSiteIndex(int s) {
		// System.err.println("setStoreSiteIndex\t"+s);
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
	

	private final NumberFormat formatter = NumberFormat.getNumberInstance();
	@Override
	public String toString(){
		formatter.setMaximumFractionDigits(3);
		StringBuffer sb = new StringBuffer();
		for (int p = 0; p < 4; p++) {
			for (int i = 0; i < getLength(); i++) {
				String freq = formatter.format(getFrequency(i, p));
				sb.append(freq).append("\t");
			}
			sb.append("\n");
		}
		return sb.toString();
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
	
	public static CategorySpectrum duplicateSpectrum(CategorySpectrum oldSpectrum) {

//		super("Spectrum");
		CategorySpectrum newSpectrum = new CategorySpectrum(oldSpectrum.getLength());
//		spectrum = new ArrayList<SpectraParameter>();
		Taxon newTaxon = oldSpectrum.getTaxon();
		newSpectrum.setTaxon(newTaxon);
		for (int i = 0; i < oldSpectrum.getLength(); i++) {
			int cat = oldSpectrum.getSpectra(i).getCategory();
//			SpectraParameter spectra = new SpectraParameter(frequencies);
			newSpectrum.getSpectra(i).setCategory(cat);
//			
//			addVariable(spectra);
//			spectrum.add(spectra);
		}
		return newSpectrum;
//		dataType = Nucleotides.INSTANCE;
//		stateCount = dataType.getStateCount();
	}
	

	public static Spectrum duplicateSpectrum(Spectrum oldSpectrum) {

		Spectrum newSpectrum = new Spectrum(oldSpectrum.getLength());

		Taxon newTaxon = oldSpectrum.getTaxon();
		newSpectrum.setTaxon(newTaxon);
		for (int i = 0; i < oldSpectrum.getLength(); i++) {
			double[] frequencies = oldSpectrum.getFrequenciesAt(i);
			newSpectrum.resetSpectra(i, frequencies);
			
		}
		return newSpectrum;
//		dataType = Nucleotides.INSTANCE;
//		stateCount = dataType.getStateCount();
	}
	
	
	
}