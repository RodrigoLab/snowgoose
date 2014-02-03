package srp.spectrum;

import java.util.Arrays;

import srp.spectrum.SpectrumAlignmentModel.Type;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

//public class Spectra implements Parameter{
public class SpectraParameter extends Parameter.Default{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -3519136356577040343L;
	public static final double[] EQUAL_FREQ = new double[]{0.25, 0.25, 0.25, 0.25};
	public static final int DIMENSION = 4;
	public static final Bounds<Double> SPECTRA_BOUNDS = new DefaultBounds(1.0, 0.0, DIMENSION);
	
	public SpectraParameter(Type type){
		this(type.getCode());
	}
	
	public SpectraParameter(int type){
		this(EQUAL_FREQ);
		if(type==1){//1 0 0 0
			double[] freq = new double[DIMENSION];
			int dim = MathUtils.nextInt(DIMENSION);
			freq[dim]=1;
			setFrequenciesQuietly(freq);
		}
		else if(type==2){//Random
			double[] freq = new double[DIMENSION];
			double sum = 0;
			for (int i = 0; i < freq.length; i++) {
				freq[i] = MathUtils.nextInt(100);
				sum += freq[i];
			}
			for (int i = 0; i < freq.length; i++) {
				freq[i] /= sum;
			}
			setFrequenciesQuietly(freq);
			
		}
		if(type==3){//0.85 0.05 0.05 0.05
			double[] freq = new double[DIMENSION];
			Arrays.fill(freq, 0.05);
			int dim = MathUtils.nextInt(DIMENSION);
			freq[dim] = 0.85;
			setFrequenciesQuietly(freq);
		}


	}
	
	
//	public SpectraParameter(){
//		this(EQUAL_FREQ);
//	}

    public SpectraParameter(double[] frequencies) {
    	super(frequencies);
    	setId("spectraSSS");
    	

        double sum = getSumOfFrequencies(frequencies);
    	if(getDimension()!=DIMENSION){
    		throw new IllegalArgumentException("Frequencies should have 4 elements, frequencies.length= "+getDimension());
    	}
        if (Math.abs(sum - 1.0) > 1e-8) {
            throw new IllegalArgumentException("Frequencies do not sum to 1, they sum to " + sum);
        }
    	
		addBounds();
		if(!isWithinBounds()){
			throw new IllegalArgumentException("Frequencies out of bounds 0 < f < 1\t"+ Arrays.toString(frequencies)); 
		}

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
    
    //Sets the value of the parameter without firing a changed event.
    protected void setFrequenciesQuietly(double[] values){
    	for (int i = 0; i < DIMENSION; i++) {
    		setParameterValueQuietly(i, values[i]);
		}
//    	fireParameterChangedEvent(i, Parameter.ChangeType.VALUE_CHANGED);
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
    	
//        double[] frequencies = new double[getFrequencyCount()];
//        for (int i = 0; i < frequencies.length; i++) {
//            frequencies[i] = getFrequency(i);
//        }
        return getParameterValues();
    }

    public double[] getCumulativeFrequencies() {
        double[] frequencies = getFrequencies();
        for (int i = 1; i < frequencies.length; i++) {
            frequencies[i] += frequencies[i - 1];
        }
        return frequencies;
    }
    
    protected void storeState() {
//    	System.err.println("storeState in Spectra");
    	super.storeValues();
	}
    protected void restoreState() {
//    	System.err.println("restoreState in Spectra");
    	super.restoreValues();
	}
    public String diagnostic(){
    	String diag = Arrays.toString(getFrequencies());
    	return diag;
    }
    
    private void addBounds(){
		addBounds(SPECTRA_BOUNDS);
	}
	
}
