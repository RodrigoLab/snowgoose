package srp.spectrum;
import java.util.ArrayList;
import java.util.Iterator;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.util.Attributable;
import dr.util.Identifiable;

//public class Spectrum implements Identifiable, Attributable{// extends AbstractModel{
public class Spectrum extends AbstractModel implements Attributable {
//TODO should be a model as well?


	private static final long serialVersionUID = -728370884996776301L;

//	private Parameter frequencyParameter;
	private ArrayList<SpectraParameter> spectrum;
	private ArrayList<SpectraParameter> storeSpectrum;
	
    protected Taxon taxon = null;
//    private StringBuffer sequenceString = null;
    protected DataType dataType = null;
    protected int stateCount;

//	public Spectra() {
//		sequenceString = new StringBuilder();
//	}

public Spectrum(int length){
			this(length, SpectraParameter.EQUAL_FREQ);
	//		spectrum = new ArrayList<Spectra>();
	//		for (int i = 0; i < length; i++) {
	//			Spectra f = new Spectra();
	//			spectrum.add(f);
	//		}
	//	    double sum = getSumOfFrequencies(frequencyParameter);
	//	
	//	    if (Math.abs(sum - 1.0) > 1e-8) {
	//	        throw new IllegalArgumentException("Frequencies do not sum to 1, they sum to " + sum);
	//	    }
	//	
	//	    this.frequencyParameter = frequencyParameter;
	//	    addVariable(frequencyParameter);
	//	    frequencyParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, frequencyParameter.getDimension()));
		    
		}

	//	public Spectra(String sequence) {
//		this();
//		setSequenceString(sequence.toUpperCase());
//	}
//
//	public Spectra(Sequence sequence) {
//		this(sequence.getTaxon(), sequence.getSequenceString());
//	}
//
//	public Spectra(Taxon taxon, String sequence) {
//		this();
//		setTaxon(taxon);
//		setSequenceString(sequence);
//	}
	public Spectrum(int length, double[] initFreq){
		
		super("Spectrum");
		spectrum = new ArrayList<SpectraParameter>();
		for (int i = 0; i < length; i++) {
			SpectraParameter spectra = new SpectraParameter(initFreq);
			addVariable(spectra);
			spectrum.add(spectra);
		}
		dataType = Nucleotides.INSTANCE;
		stateCount = dataType.getStateCount();
	}
	
	public Spectrum(Sequence sequence) {
			this(sequence.getLength());
	//		setDataType(sequence.getDataType());
	//		super("Spectrum");
			dataType = Nucleotides.INSTANCE;
			stateCount = dataType.getStateCount();
			for (int site = 0; site < sequence.getLength(); site++) {
				char c = sequence.getChar(site);
				int state = dataType.getState(c);
				if (state < stateCount) {
					double[] newFreq = new double[stateCount];
					newFreq[state]=1;
					setFrequencies(site, newFreq);
				}
				else{
					setFrequencyAt(site, 0, 0.25);
					setFrequencyAt(site, 1, 0.25);
					setFrequencyAt(site, 2, 0.25);
					setFrequencyAt(site, 3, 0.25);
				}
			}
		
		}

	//	/**
//	 * @param frequencies the frequencies
//	 * @return return the sum of frequencies
//	 */
//	private double getSumOfFrequencies(Parameter frequencies) {
//	    double total = 0.0;
//	    for (int i = 0; i < frequencies.getDimension(); i++) {
//	        total += frequencies.getParameterValue(i);
//	    }
//	    return total;
//	}
	public int getLength(){
		return spectrum.size();
	}
	public void setFrequencyAt(int site, int state, double value) {
		spectrum.get(site).setFrequency(state, value);
		fireModelChanged(this);//TODO check?
	}
	public void setFrequencies(int site, double[] values){
		spectrum.get(site).setFrequencies(values);
	}
	public double getFrequency(int site, int state) {
	    return spectrum.get(site).getFrequency(state);
	}
	
	public double[] getFrequencies(int site) {
		return getSpecturm(site).getFrequencies();
	}

	public SpectraParameter getSpecturm(int site) {
	    return spectrum.get(site);
	}
	
	@Override
	public String toString(){
		StringBuffer sb = new StringBuffer();
		for (int p = 0; p < 4; p++) {
			for (int i = 0; i < spectrum.size(); i++) {
				sb.append(getFrequency(i, p)).append("\t");
			}
			sb.append("\n");
		}
		return sb.toString();
	}
//	public int getFrequencyCount() {
//	    return frequencyParameter.getDimension();
//	}
//	
//	public Parameter getFrequencyParameter() {
//	    return frequencyParameter;
//	}
	
//	public double[] getFrequencies() {
//	    double[] frequencies = new double[getFrequencyCount()];
//	    for (int i = 0; i < frequencies.length; i++) {
//	        frequencies[i] = getFrequency(i);
//	    }
//	    return frequencies;
//	}
	
//	public double[] getCumulativeFrequencies() {
//	    double[] frequencies = getFrequencies();
//	    for (int i = 1; i < frequencies.length; i++) {
//	        frequencies[i] += frequencies[i - 1];
//	    }
//	    return frequencies;
//	}

	
	
	
	
	
	
	
	
//	public void setCharAt(int index, int newChar) {
//		sequenceString.setCharAt(index, (char) newChar);
//	}
//	
//	public void setCharAt(int index, char newChar) {
//		sequenceString.setCharAt(index, newChar);
//	}
//
//	public char replaceCharAt(int index, int newChar){
//		char oldChar = getChar(index);
//		setCharAt(index, (char) newChar);
//		return oldChar;
//	}
//	
	// **************************************
	// OVERRIDE ALL (almost all) methods
	// Do NOT call setState()!!
	// ************************************
	



//    /**
//     * @return the length of the sequences.
//     */
//    @Override
//	public int getLength() {
//        return sequenceString.length();
//    }
//
//    /**
//     * @return a String containing the sequences.
//     */
//    @Override
//	public String getSequenceString() {
//        return sequenceString.toString();
//    }
//
//    /**
//     * @return a char containing the state at index.
//     */
//    @Override
//	public char getChar(int index) {
//        return sequenceString.charAt(index);
//    }
//
//    /**
//     * @return the state at site index.
//     */
//    @Override
//	public int getState(int index) {
//        return dataType.getState(sequenceString.charAt(index));
//    }
//
//    /**
//     */
////		public void setState(int index, int state) {
////
////	        sequenceString.setCharAt(index, dataType.getChar(state));
////	    }
//
//    /**
//     * Characters are copied from the sequences into the destination character array dst.
//     */
//    @Override
//	public void getChars(int srcBegin, int srcEnd, char[] dst, int dstBegin) {
//        sequenceString.getChars(srcBegin, srcEnd, dst, dstBegin);
//    }
//
//    /**
//     * Set the DataType of the sequences.
//     */
//    @Override
//	public DataType guessDataType() {
//        return DataType.guessDataType(sequenceString.toString());
//    }
//
//    /**
//     * Set the sequences using a string.
//     */
//    @Override
//	public void setSequenceString(String sequence) {
//        sequenceString.setLength(0);
//        sequenceString.append(sequence);
//    }
//
//    /**
//     * Append a string to the sequences.
//     */
//    @Override
//	public void appendSequenceString(String sequence) {
//        sequenceString.append(sequence);
//    }
//
//    /**
//     * Insert a string into the sequences.
//     */
//    @Override
//	public void insertSequenceString(int offset, String sequence) {
//        sequenceString.insert(offset, sequence);
//    }

    
    public void setDataType(DataType dataType) {
        this.dataType = dataType;
    }

	public DataType getDataType() {
        return dataType;
    }
	


/////////////////////////

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
    // **************************************************************
    // AbstractModel IMPLEMENTATION
    // **************************************************************


	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		System.err.println("Call handleModelChangedEvent");
		
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		System.err.println("Call handleVariableChangedEvent Specturm");
		
	}

	@Override
	public void storeState() {
		//TODO don't need to copy everything
		for (int i = 0; i < spectrum.size(); i++) {
			storeSpectrum.set(i, spectrum.get(i));
		}
	}

	@Override
	public void restoreState(){
		//TODO don't need to copy everything
		for (int i = 0; i < spectrum.size(); i++) {
			spectrum.set(i, storeSpectrum.get(i));
		}

	}

	@Override
	protected void acceptState() {
		// Do nothing
		
	}

    	
}
