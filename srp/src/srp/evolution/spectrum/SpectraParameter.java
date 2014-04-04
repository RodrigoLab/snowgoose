package srp.evolution.spectrum;

import java.util.Arrays;

import srp.dr.evolution.datatype.ShortReads;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

public class SpectraParameter extends AbstractSpectra{
	
	private static final long serialVersionUID = 3708765314559863330L;

	public static final int TOTAL_STATE_COUNT = ShortReads.INSTANCE.getAmbiguousStateCount();
	public static final int DIMENSION = 4;
	public static final double[] EQUAL_FREQ = new double[]{0.25, 0.25, 0.25, 0.25};
	public static final Bounds<Double> SPECTRA_BOUNDS = new DefaultBounds(1.0, 0.0, DIMENSION);

	public static final double MAX_FREQ = 0.99;
	public static final double MIN_FREQ = 0.001;//MrBayes 0.0001

	public static final double INIT_LARGE = 0.97;
	public static final double INIT_SMALL = 0.01;

	
    protected double[] values;
    protected double[] storedValues;
    
	protected double[] storedstateLikelihood;
	protected double[] stateLikelihood;

	
	
	public SpectraParameter(SpectraType type){
		this(EQUAL_FREQ);
		double[] freq = new double[DIMENSION];
		
		switch (type) {
		case DOMINANT:
			Arrays.fill(freq, INIT_SMALL);
			int dim = MathUtils.nextInt(DIMENSION);
			freq[dim] = INIT_LARGE;
			break;
		case RANDOM:
			double sum = 0;
			for (int i = 0; i < freq.length; i++) {
				freq[i] = 0.01 + MathUtils.nextDouble();// +MathUtils.nextInt(1000);
				sum += freq[i];
			}
			for (int i = 0; i < freq.length; i++) {
				freq[i] /= sum;
			}
			break;
		case CATEGORY:
			Arrays.fill(freq, 0.1);
			freq[MathUtils.nextInt(DIMENSION)] = 0.7;
		case EQUAL:
			Arrays.fill(freq, 0.25);
			break;
		default:
			throw new IllegalArgumentException("Invalid type: "+type);
		}
		setFrequenciesQuietly(freq);

	}
	

	
    public SpectraParameter(double[] frequencies) {

    	setId("spectra");
    	
    	dimension = DIMENSION;
    	if(frequencies.length!=DIMENSION){
    		throw new IllegalArgumentException("Frequencies should have 4 elements, frequencies.length= "+frequencies.length);
    	}
        addBounds(SPECTRA_BOUNDS);
        this.values = new double[DIMENSION];
        this.storedValues = new double[DIMENSION];
        
        this.stateLikelihood = new double[TOTAL_STATE_COUNT];
        this.storedstateLikelihood = new double[TOTAL_STATE_COUNT];
        
        System.arraycopy(frequencies, 0, values, 0, DIMENSION);
        checkSpectra();

    }

    private void setFrequenciesQuietly(double[] values){
		for (int i = 0; i < DIMENSION; i++) {
			setFrequency(i, values[i]);
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
        values[i] = value;
        fireParameterChangedEvent(i, Parameter.ChangeType.VALUE_CHANGED);

    }
    
    public double getFrequency(int i) {
    	return values[i];
    }

    public double[] getFrequencies() {
    	
    	double[] copyOfValues = new double[DIMENSION];
        System.arraycopy(values, 0, copyOfValues, 0, DIMENSION);
        return copyOfValues;
    }

	@Override
    public boolean isWithinBounds() {
        Bounds<Double> bounds = getBounds();
        for (int i = 0; i < getDimension(); i++) {
            final double value = getFrequency(i);
            if (value < bounds.getLowerLimit(i) || value > bounds.getUpperLimit(i)) {
                return false;
            }
        }
        return true;
    }
    public void setStateLikelihood(double[] likelihood){
    	System.arraycopy(likelihood, 0, stateLikelihood, 0, DIMENSION);
    }
    
    public void setStateLikelihood(double[] allStateLikelihood, int offset){
    	System.arraycopy(allStateLikelihood, offset, stateLikelihood, 0, DIMENSION);
    }
	public void setStateLikelihood(int state, double ln) {
		stateLikelihood[state] = ln;
	}
	
	public double[] getStateLikelihood() {
		return stateLikelihood;
	}
    public double[] getStoredStateLikelihood() {
		return storedstateLikelihood;
	}
    
	public void storeState() {
//    	System.err.println("storeState in Spectra");
    	storeValues();
	}
    public void restoreState() {
//    	System.err.println("restoreState in Spectra");
    	restoreValues();
	}


    public void checkSpectra(){
//		double sum = 0;
//		for (int j = 0; j < DIMENSION; j++) {
//			double f = getFrequency(j);
//			sum += f;
//			if(f<0 || f>1){
//				System.err.println(j +"\t"+ f +"\t"+ Arrays.toString(sp.getFrequencies()));
//			}
//			
//		}
//		if(sum>1.01 || sum<0.99){
//			System.err.println(Arrays.toString(sp.getFrequencies()));
//		}
		
		 double sum = getSumOfFrequencies(values);
	    	
	        if (Math.abs(sum - 1.0) > 1e-8) {
	            throw new IllegalArgumentException("Frequencies do not sum to 1, they sum to " + sum);
	        }

			if(!isWithinBounds()){
				throw new IllegalArgumentException("Frequencies out of bounds 0 < f < 1\t"+ Arrays.toString(values)); 
			}
	}
    
    
	@Deprecated
	public double getStoredFrequency(int state) {
		return storedValues[state];
	}


	
	public enum SpectraType{
		DOMINANT,
		RANDOM,
		@Deprecated EQUAL,
		@Deprecated CATEGORY,

	}



    @Override
	protected final void storeValues() {
        System.arraycopy(values, 0, storedValues, 0, DIMENSION);
        System.arraycopy(stateLikelihood, 0, storedstateLikelihood, 0, DIMENSION);
    }

    @Override
	protected final void restoreValues() {
        //swap the arrays
        double[] temp = storedValues;
        storedValues = values;
        values = temp;
        
        temp = storedstateLikelihood;
    	storedstateLikelihood = stateLikelihood;
    	stateLikelihood = temp;

    }

    /**
     * Nothing to do
     */
    @Override
	protected final void acceptValues() {
    }




}
//
//
