package srp.spectrum;

import java.util.Arrays;

import srp.dr.evolution.datatype.ShortReads;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

public class SpectraParameter extends AbstractSpectra{
	
	private static final long serialVersionUID = 3708765314559863330L;

	public static final int TOTAL_STATE_COUNT = ShortReads.INSTANCE.getAmbiguousStateCount();
	public static final double MAX_FREQ = 0.99;
	public static final double MIN_FREQ = 0.001;//MrBayes 0.0001
	
	private double[] storedstateLikelihood;
	private double[] stateLikelihood;
	
	public SpectraParameter(SpectraType type){
		this(EQUAL_FREQ);
		double[] freq = new double[DIMENSION];
		switch (type) {
		case ZERO_ONE:
			Arrays.fill(freq, 0.01);
			int dim = MathUtils.nextInt(DIMENSION);
			freq[dim]=0.97;
			setFrequenciesQuietly(freq);
			break;
		case RANDOM:
			double sum = 0;
			for (int i = 0; i < freq.length; i++) {
				freq[i] = 0.01+MathUtils.nextDouble();//+MathUtils.nextInt(1000);
				sum += freq[i];
			}
			for (int i = 0; i < freq.length; i++) {
				freq[i] /= sum;
			}
			setFrequenciesQuietly(freq);
			break;
		case CATEGORY:
			Arrays.fill(freq, 0.1);
			freq[MathUtils.nextInt(DIMENSION)] = 0.7;
			setFrequenciesQuietly(freq);
		default:
			break;
		}
	}
	

	
    public SpectraParameter(double[] frequencies) {
    	
    	if(frequencies.length!=DIMENSION){
    		throw new IllegalArgumentException("Frequencies should have 4 elements, frequencies.length= "+frequencies.length);
    	}
        addBounds(SPECTRA_BOUNDS);
        this.values = new double[DIMENSION];
        this.storedValues = new double[DIMENSION];
        
        this.stateLikelihood = new double[TOTAL_STATE_COUNT];
        this.storedstateLikelihood = new double[TOTAL_STATE_COUNT];
        
        System.arraycopy(frequencies, 0, values, 0, DIMENSION);

    	setId("spectra");
    	

        double sum = getSumOfFrequencies(frequencies);
    	
        if (Math.abs(sum - 1.0) > 1e-8) {
            throw new IllegalArgumentException("Frequencies do not sum to 1, they sum to " + sum);
        }
    	

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
    
    protected void setFrequenciesQuietly(double[] values){
    	for (int i = 0; i < DIMENSION; i++) {
    		setParameterValueQuietly(i, values[i]);
		}
    }
    
    public double getFrequency(int i) {
        return getParameterValue(i);
    }

    public double[] getFrequencies() {
        return getParameterValues();
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
	@Deprecated
	public double[] getStateLikelihood() {
		return stateLikelihood;
	}
    public double[] getStoredStateLikelihood() {
		return storedstateLikelihood;
	}

	public void storeState() {
//    	System.err.println("storeState in Spectra");
    	super.storeValues();
    	
//        System.arraycopy(stateLikelihood, 0, storedstateLikelihood, 0, stateLikelihood.length);
        System.arraycopy(stateLikelihood, 0, storedstateLikelihood, 0, DIMENSION);
//        System.out.println("STORE:\t"+Arrays.toString(stateLikelihood) +"\t"+ Arrays.toString(storedstateLikelihood));
	}
    public void restoreState() {
//    	System.err.println("restoreState in Spectra");
    	super.restoreValues();
    	double[] temp = storedstateLikelihood;
    	storedstateLikelihood = stateLikelihood;
    	stateLikelihood = temp;
//    	System.arraycopy(storedstateLikelihood, 0, stateLikelihood, 0, DIMENSION);
	}


	@Deprecated
	public double getStoredFrequency(int state) {
		return storedValues[state];
	}


	
	public enum SpectraType{
		@Deprecated EQUAL,
		@Deprecated ZERO_ONE,
		RANDOM, 
		@Deprecated CATEGORY;

	}




}

abstract class AbstractSpectra extends Parameter.Abstract {		
	

	public static final double[] EQUAL_FREQ = new double[]{0.25, 0.25, 0.25, 0.25};
	public static final int DIMENSION = 4;
	public static final Bounds<Double> SPECTRA_BOUNDS = new DefaultBounds(1.0, 0.0, DIMENSION);

	private static final long serialVersionUID = 8795190465907061489L;
	
    protected double[] values;
    protected double[] storedValues;

    private Bounds<Double> bounds = null;
    


	@Override
	public void addBounds(Bounds<Double> boundary) {
        if (bounds == null) {
            bounds = boundary;
        } else {
            throw new IllegalArgumentException("Should not call addBounds twice");
            // can't change dimension after bounds are added!
        }

    }

    //********************************************************************
    // GETTERS
    //********************************************************************

    @Override
	public final int getDimension() {
        return DIMENSION;
    }

    @Override
	public final int getSize() {
        return DIMENSION;
    }

    @Override
	public final double getParameterValue(int i) {
        return values[i];
    }

    /**
     * Defensively returns copy of parameter array.
     *
     * @return a copy of the parameter values
     */
    @Override
	public final double[] getParameterValues() {

        double[] copyOfValues = new double[DIMENSION];
        System.arraycopy(values, 0, copyOfValues, 0, DIMENSION);
        return copyOfValues;
    }

    @Override
	public Bounds<Double> getBounds() {
        if (bounds == null) {
            throw new NullPointerException(getParameterName() + " parameter: Bounds not set");
        }
        return bounds;
    }

    @Override
	public String getParameterName() {
        return getId();
    }

    //********************************************************************
    // SETTERS
    //********************************************************************



    @Override
	public void setParameterValue(int i, double val) {
        values[i] = val;
        fireParameterChangedEvent(i, Parameter.ChangeType.VALUE_CHANGED);
    }

    /**
     * Sets the value of the parameter without firing a changed event.
     *
     * @param dim   the index of the parameter dimension
     * @param value the value to set
     */
    @Override
	public void setParameterValueQuietly(int dim, double value) {
        values[dim] = value;
    }


    /**
     * Sets the values of the parameter and notify that all values of the parameter have changed.
     *
     * @param i   index of the value
     * @param val to value to set
     */
    @Override
	public void setParameterValueNotifyChangedAll(int i, double val) {
        values[i] = val;
        fireParameterChangedEvent(i, Parameter.ChangeType.ALL_VALUES_CHANGED);
    }

    @Override
	protected final void storeValues() {
        System.arraycopy(values, 0, storedValues, 0, DIMENSION);
    }

    @Override
	protected final void restoreValues() {
        //swap the arrays
        double[] temp = storedValues;
        storedValues = values;
        values = temp;

    }

    /**
     * Nothing to do
     */
    @Override
	protected final void acceptValues() {
    }

    @Override
	protected final void adoptValues(Parameter source) {
    	throw new IllegalArgumentException("Can not adoptValues");
    }
    
    @Override
	public void setDimension(int dim) {
    	throw new IllegalArgumentException("Can not setDimension");
    }

    @Override
	public void addDimension(int index, double value) {
    	throw new IllegalArgumentException("Can not addDimension");
    }
    @Override
	public double removeDimension(int index) {
    	throw new IllegalArgumentException("Can not removeDimension");
    }


}
//////////////////////////////
//
//package srp.spectrum;
//
//import java.util.Arrays;
//
//import dr.inference.model.Bounds;
//import dr.inference.model.Parameter;
//import dr.math.MathUtils;
//
////public class Spectra implements Parameter{
//public class SpectraParameter extends Parameter.Default{
//	
//
//	private static final long serialVersionUID = -3519136356577040343L;
//	public static final double[] EQUAL_FREQ = new double[]{0.25, 0.25, 0.25, 0.25};
//	public static final int DIMENSION = 4;
//	public static final Bounds<Double> SPECTRA_BOUNDS = new DefaultBounds(1.0, 0.0, DIMENSION);
//	public static final double MAX_FREQ = 0.99;
//	public static final double MIN_FREQ = 0.001;//MrBayes 0.0001
//	
//	
//	public SpectraParameter(SpectraType type){
//		this(EQUAL_FREQ);
//		double[] freq = new double[DIMENSION];
//		switch (type) {
//		case ZERO_ONE:
//			int dim = MathUtils.nextInt(DIMENSION);
//			freq[dim]=1;
//			setFrequenciesQuietly(freq);
//			break;
//		case RANDOM:
//			double sum = 0;
//			for (int i = 0; i < freq.length; i++) {
//				freq[i] = 0.01+MathUtils.nextDouble();//+MathUtils.nextInt(1000);
//				sum += freq[i];
//			}
//			for (int i = 0; i < freq.length; i++) {
//				freq[i] /= sum;
//				if(freq[i]<=0){
//					System.err.println(Arrays.toString(freq) +"\t"+ sum);
//				}
//			}
//			
//			setFrequenciesQuietly(freq);
//			break;
//		case CATEGORY:
//			Arrays.fill(freq, 0.05);
//			freq[MathUtils.nextInt(DIMENSION)] = 0.85;
//			setFrequenciesQuietly(freq);
//		default:
//			break;
//		}
//	}
//	
//	@Deprecated
//	public SpectraParameter(int type){
//		this(EQUAL_FREQ);
//		
//
//
//	}
//	
//	
//    public SpectraParameter(double[] frequencies) {
//    	super(frequencies);
//    	setId("spectraSSS");
//    	
//
//        double sum = getSumOfFrequencies(frequencies);
//    	if(getDimension()!=DIMENSION){
//    		throw new IllegalArgumentException("Frequencies should have 4 elements, frequencies.length= "+getDimension());
//    	}
//        if (Math.abs(sum - 1.0) > 1e-8) {
//            throw new IllegalArgumentException("Frequencies do not sum to 1, they sum to " + sum);
//        }
//    	
//		addBounds();
//		if(!isWithinBounds()){
//			throw new IllegalArgumentException("Frequencies out of bounds 0 < f < 1\t"+ Arrays.toString(frequencies)); 
//		}
//
//    }
//
//    public void setFrequency(int i, double value) {
//    	setParameterValue(i, value);
//    }
//    
//    //Sets the value of the parameter without firing a changed event.
//    protected void setFrequenciesQuietly(double[] values){
//    	for (int i = 0; i < DIMENSION; i++) {
//    		setParameterValueQuietly(i, values[i]);
//		}
////    	fireParameterChangedEvent(i, Parameter.ChangeType.VALUE_CHANGED);
////    	System.arraycopy(values, 0, this.values, 0, values.length);
//    }
//    
//    public double getFrequency(int i) {
//        return getParameterValue(i);
//    }
//
//    public int getFrequencyCount() {
//        return getDimension();
//    }
//
////    public Parameter getFrequencyParameter() {
////        return spectra;
////    }
//
//    public double[] getFrequencies() {
//    	
////        double[] frequencies = new double[getFrequencyCount()];
////        for (int i = 0; i < frequencies.length; i++) {
////            frequencies[i] = getFrequency(i);
////        }
//        return getParameterValues();
//    }
//
//
//    
//    public void storeState() {
////    	System.err.println("storeState in Spectra");
//    	super.storeValues();
//    	
//        System.arraycopy(getParameterValues(), 0, storedFrequency, 0, storedFrequency.length);
//    
//	}
//    protected void restoreState() {
////    	System.err.println("restoreState in Spectra");
//    	super.restoreValues();
//	}
//
//    private void addBounds(){
//		addBounds(SPECTRA_BOUNDS);
//	}
//
//	public double getStoredFrequency(int state) {
////		storedValues
//		
//		return storedFrequency[state];
//	}
//	private double[] storedFrequency = new double[4];
//
//	
//	private static double getSumOfFrequencies(double[] frequencies) {
//	    double total = 0.0;
//	    for (int i = 0; i < frequencies.length; i++) {
//	        total += frequencies[i];
//	    }
//	    return total;
//	}
//
//
//	public enum SpectraType{
//		EQUAL,
//		ZERO_ONE,
//		RANDOM, 
//		CATEGORY;
//
//	}
//
//}
//
//
