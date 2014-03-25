package srp.spectrum;
import java.text.NumberFormat;
import java.util.ArrayList;

import srp.spectrum.CategorySpectraParameter.CategoryType;
import dr.evolution.util.Taxon;

//public class Spectrum implements Identifiable, Attributable{// extends AbstractModel{
public class CategorySpectrum extends AbstractSpectrum {

	private static final long serialVersionUID = -728370884996776301L;

	ArrayList<CategorySpectraParameter> spectrum;

	private CategorySpectraParameter[] spectrumArray;
    
//	public Spectra() {
//		sequenceString = new StringBuilder();
//	}

	public CategorySpectrum(int length) {
		this(length, CategoryType.SINGLE);
	}


//	public CategorySpectrum(double[][] freqs) {
//		this(freqs[0].length);
//		for (int s = 0; s < spectrum.size(); s++) {
//			CategorySpectraParameter spectra = spectrum.get(s);
//			for (int f = 0; f < 4; f++) {
//				spectra.setFrequency(f, freqs[f][s]);
//			}
////			double[] freq = new double[]{freqs[0][s], freqs[1][s], freqs[2][s], freqs[3][s]}; 
////			spectrum.resetFrequencies(s, freq);
//		}
//	}

	/**
	 * @param length
	 * @param type
	 */
	public CategorySpectrum(int length, CategoryType type){
		
		super("CategorySpectrum");
		spectrumArray = new CategorySpectraParameter[length];
		spectrum = new ArrayList<CategorySpectraParameter>();
		for (int i = 0; i < length; i++) {
			CategorySpectraParameter spectra = new CategorySpectraParameter(type);
			addVariable(spectra);
			spectrum.add(spectra);
			spectrumArray[i] = spectra;
		}
		
	}

	
//	public CategorySpectrum(Sequence sequence) {
//		this(sequence.getLength());
//		
//		for (int i = 0; i < sequence.getLength(); i++) {
//			SpectraParameter spectra = getSpectra(i);
//			char c = sequence.getChar(i);
//			int state = DATA_TYPE.getState(c);
//			for (int j = 0; j < STATE_COUNT; j++) {
//				if(j==state){
//					spectra.setFrequency(j, 1);
//				}
//				else{
//					spectra.setFrequency(j, 0);
//				}
//			}
//		}
//		
//		}

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


	
	public int getLength() {
		return spectrum.size();
	}


//	public void resetFrequencies(int site, double[] values){
//		getSpectra(site).setFrequenciesQuietly(values);
//	}
	
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

	public CategorySpectraParameter getSpectra(int site) {
	    return spectrum.get(site);
	}

	public CategorySpectraParameter getSpectraArray(int site) {
	    return spectrumArray[site];
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

    

/////////////////////////


    	
}
