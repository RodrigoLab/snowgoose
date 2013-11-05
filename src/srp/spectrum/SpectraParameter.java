package srp.spectrum;

import java.util.Arrays;

import dr.inference.model.Parameter;

//public class Spectra implements Parameter{
public class SpectraParameter extends Parameter.Default{
	
	public static final double[] EQUAL_FREQ = new double[]{0.25, 0.25, 0.25, 0.25};
	
	public SpectraParameter(){
		this(EQUAL_FREQ);
	}

    public SpectraParameter(double[] frequencies) {
    	super(frequencies);
    	

        double sum = getSumOfFrequencies(frequencies);
    	if(getDimension()!=4){
    		throw new IllegalArgumentException("Frequencies should have 4 elements, frequencies.length= "+getDimension());
    	}
        if (Math.abs(sum - 1.0) > 1e-8) {
            throw new IllegalArgumentException("Frequencies do not sum to 1, they sum to " + sum);
        }
    	
		addBounds(new DefaultBounds(1.0, 0.0, getDimension()));
		if(!isWithinBounds()){
			throw new IllegalArgumentException("Frequencies out of bounds 0 < f < 1\t"+ Arrays.toString(frequencies)); 
		}
//        this.frequencyParameter = frequencyParameter;
//        addVariable(frequencyParameter);
//        frequencyParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, frequencyParameter.getDimension()));
//        this.dataType = dataType;
    }

    private static double getSumOfFrequencies(double[] frequencies) {
        double total = 0.0;
        for (int i = 0; i < frequencies.length; i++) {
            total += frequencies[i];
        }
        return total;
    }

    public void setFrequency(int i, double value) {
    	setParameterValue(i, value);
    }
    
    public void setFrequencies(double[] values){
    	for (int i = 0; i < values.length; i++) {
    		setParameterValueQuietly(i, values[i]);
		}
//    	System.arraycopy(values, 0, this.values, 0, values.length);
    }
    
    public double getFrequency(int i) {
        return getParameterValue(i);
    }

    public int getFrequencyCount() {
        return getDimension();
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
