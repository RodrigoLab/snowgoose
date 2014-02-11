package srp.spectrum;
import java.text.NumberFormat;
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

//public class Spectrum implements Identifiable, Attributable{// extends AbstractModel{
public class Spectrum extends AbstractModel implements Attributable {



	private static final long serialVersionUID = -728370884996776301L;

//	private Parameter frequencyParameter;
	private ArrayList<SpectraParameter> spectrum;
	private ArrayList<SpectraParameter> storeSpectrum;
	
    protected Taxon taxon = null;
//    private StringBuffer sequenceString = null;
    protected DataType dataType = null;
    protected int stateCount;
    private int storeSiteIndex;
    
//	public Spectra() {
//		sequenceString = new StringBuilder();
//	}

	public Spectrum(int length) {
		this(length, 0);
	}


	public Spectrum(double[][] freqs) {
		this(freqs[0].length);
		for (int s = 0; s < spectrum.size(); s++) {
			SpectraParameter spectra = spectrum.get(s);
			for (int f = 0; f < 4; f++) {
				spectra.setFrequency(f, freqs[f][s]);
			}
//			double[] freq = new double[]{freqs[0][s], freqs[1][s], freqs[2][s], freqs[3][s]}; 
//			spectrum.resetFrequencies(s, freq);
		}
	}

	/**
	 * type 0=Equal. 
	 * 		1=[1 0 0 0].
	 * 		2=[Random]. 
	 * @param length
	 * @param type
	 */
	public Spectrum(int length, int type){
		
		super("Spectrum");
		spectrum = new ArrayList<SpectraParameter>();
		for (int i = 0; i < length; i++) {
			SpectraParameter spectra = new SpectraParameter(type);
			addVariable(spectra);
			spectrum.add(spectra);
		}
		dataType = Nucleotides.INSTANCE;
		stateCount = dataType.getStateCount();
	}
	
//	public static SpectrumFactory(Sequence sequence)
	public Spectrum(Sequence sequence) {
		//TODO change to static factory method
		super("Spectrum");
		
		dataType = Nucleotides.INSTANCE;
		stateCount = dataType.getStateCount();
		
		spectrum = new ArrayList<SpectraParameter>();
		SpectraParameter spectra;
		for (int i = 0; i < sequence.getLength(); i++) {
			char c = sequence.getChar(i);
			int state = dataType.getState(c);
			if (state < stateCount) {
				double[] newFreq = new double[stateCount];
				newFreq[state]=1;
				spectra = new SpectraParameter(newFreq);
//				setFrequencies(site, newFreq);
			}
			else{
				System.out.println(c +"\t"+ state +"\tat:"+ i);
				spectra = new SpectraParameter(SpectraParameter.EQUAL_FREQ);
				spectra.setFrequenciesQuietly(new double[]{1,1,1,1});
			}
			
//			SpectraParameter spectra = new SpectraParameter(initFreq);
			addVariable(spectra);
			spectrum.add(spectra);
		}
//		dataType = Nucleotides.INSTANCE;
//		stateCount = dataType.getStateCount();
		
		
		
		
		
		
		
//			this(sequence.getLength());
	//		setDataType(sequence.getDataType());
	//		super("Spectrum");
//			dataType = Nucleotides.INSTANCE;
//			stateCount = dataType.getStateCount();
//			for (int site = 0; site < sequence.getLength(); site++) {
//				char c = sequence.getChar(site);
//				int state = dataType.getState(c);
//				if (state < stateCount) {
//					double[] newFreq = new double[stateCount];
//					newFreq[state]=1;
//					setFrequencies(site, newFreq);
//				}
//				else{
//					setFrequencyAt(site, 0, 0.25);
//					setFrequencyAt(site, 1, 0.25);
//					setFrequencyAt(site, 2, 0.25);
//					setFrequencyAt(site, 3, 0.25);
//				}
//			}
		
		}

	public static Spectrum duplicateSpectrum(Spectrum oldSpectrum) {

//		super("Spectrum");
		Spectrum newSpectrum = new Spectrum(oldSpectrum.getLength());
//		spectrum = new ArrayList<SpectraParameter>();
		Taxon newTaxon = oldSpectrum.getTaxon();
		newSpectrum.setTaxon(newTaxon);
		for (int i = 0; i < oldSpectrum.getLength(); i++) {
			double[] frequencies = oldSpectrum.getFrequenciesAt(i);
			newSpectrum.resetFrequencies(i, frequencies);
//			SpectraParameter spectra = new SpectraParameter(frequencies);
//			addVariable(spectra);
//			spectrum.add(spectra);
		}
		return newSpectrum;
//		dataType = Nucleotides.INSTANCE;
//		stateCount = dataType.getStateCount();
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
	
	@Deprecated
	public void setFrequencyAt(int site, int state, double value) {
		getSpectra(site).setFrequency(state, value);
//		fireModelChanged(this);//TODO check?
	}
	public void resetFrequencies(int site, double[] values){
		getSpectra(site).setFrequenciesQuietly(values);
	}
	
	public double getFrequency(int site, int state) {
	    return getSpectra(site).getFrequency(state);
	}
	/**
	 * Create a new double[], might be slower then getFrequency(int site, int state) 
	 * bench mark required
	 * @param site
	 * @return frequencies[]
	 */
	public double[] getFrequenciesAt(int site) {
		return getSpectra(site).getFrequencies();
	}

	public SpectraParameter getSpectra(int site) {
	    return spectrum.get(site);
	}
	
	private final NumberFormat formatter = NumberFormat.getNumberInstance();
	@Override
	public String toString(){
		formatter.setMaximumFractionDigits(3);
		StringBuffer sb = new StringBuffer();
		for (int p = 0; p < 4; p++) {
			for (int i = 0; i < spectrum.size(); i++) {
				String freq = formatter.format(getFrequency(i, p));
				sb.append(freq).append("\t");
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
//		System.err.println("Call handleVariableChangedEvent Specturm\t"+variable.getVariableName() +"\t"+ index);
		//TODO implement when secptra changed
	}

	
	public void setStoreSiteIndex(int s){
//		System.err.println("setStoreSiteIndex\t"+s);
		storeSiteIndex = s;
	}
	@Override
	protected void storeState() {
//		System.err.println("storeState Specturm: "+storeSiteIndex);

		getSpectra(storeSiteIndex).storeState();
//		for (int i = 0; i < spectrum.size(); i++) {
//			storeSpectrum.set(i, spectrum.get(i));
		
//		}
	}

	@Override
	protected void restoreState(){
//		System.err.println("restoreState Spectrum: "+storeSiteIndex);
		getSpectra(storeSiteIndex).restoreState();

//		for (int i = 0; i < spectrum.size(); i++) {
//			spectrum.set(i, storeSpectrum.get(i));
//		}

	}

	@Override
	protected void acceptState() {
		// Do nothing
		
	}

    	
}
