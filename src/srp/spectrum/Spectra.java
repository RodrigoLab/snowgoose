package srp.spectrum;

import java.util.ArrayList;

import dr.evolution.datatype.DataType;
import dr.evomodelxml.substmodel.FrequencyModelParser;
import dr.inference.model.Parameter;

//public class Spectra implements Parameter{
@Deprecated
public class Spectra {
	
	public static final double[] EQUAL_FREQ = new double[]{0.25, 0.25, 0.25, 0.25};
	
	private Parameter spectra;
	
	public Spectra(){
		
		spectra = new Parameter.Default(EQUAL_FREQ);
		spectra.addBounds(new Parameter.DefaultBounds(1.0, 0.0, spectra.getDimension()));
	}
	

    public Spectra(double[] frequencies) {

//        super(FrequencyModelParser.FREQUENCY_MODEL);
    	if(frequencies.length != 4){
    		throw new IllegalArgumentException("Frequencies should have 4 elements");
    	}
        double sum = getSumOfFrequencies(frequencies);
        if (Math.abs(sum - 1.0) > 1e-8) {
            throw new IllegalArgumentException("Frequencies do not sum to 1, they sum to " + sum);
        }

        spectra = new Parameter.Default(frequencies);
		spectra.addBounds(new Parameter.DefaultBounds(1.0, 0.0, spectra.getDimension()));
		
//        this.frequencyParameter = frequencyParameter;
//        addVariable(frequencyParameter);
//        frequencyParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, frequencyParameter.getDimension()));
//        this.dataType = dataType;
    }

    /**
     * @param frequencies the frequencies
     * @return return the sum of frequencies
     */
    private double getSumOfFrequencies(Parameter frequencies) {
        double total = 0.0;
        for (int i = 0; i < frequencies.getDimension(); i++) {
            total += frequencies.getParameterValue(i);
        }
        return total;
    }
    
    private double getSumOfFrequencies(double[] frequencies) {
        double total = 0.0;
        for (int i = 0; i < frequencies.length; i++) {
            total += frequencies[i];
        }
        return total;
    }

    public void setFrequency(int i, double value) {
    	spectra.setParameterValue(i, value);
    }

    public double getFrequency(int i) {
        return spectra.getParameterValue(i);
    }

    public int getFrequencyCount() {
        return spectra.getDimension();
    }

//    public Parameter getFrequencyParameter() {
//        return spectra;
//    }

    public double[] getFrequencies() {
        double[] frequencies = new double[getFrequencyCount()];
        for (int i = 0; i < frequencies.length; i++) {
            frequencies[i] = getFrequency(i);
        }
        return frequencies;
    }

    public double[] getCumulativeFrequencies() {
        double[] frequencies = getFrequencies();
        for (int i = 1; i < frequencies.length; i++) {
            frequencies[i] += frequencies[i - 1];
        }
        return frequencies;
    }

}
