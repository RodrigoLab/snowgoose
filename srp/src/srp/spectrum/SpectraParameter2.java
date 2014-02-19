package srp.spectrum;

import java.util.Arrays;

import srp.spectrum.SpectrumAlignmentModel.SpectrumType;
import dr.inference.model.Bounds;
import dr.inference.model.IntersectionBounds;
import dr.inference.model.Parameter;
import dr.inference.model.Parameter.DefaultBounds;
import dr.math.MathUtils;

//public class Spectra implements Parameter{
public class SpectraParameter2 extends AbstractSpectra{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -3519136356577040343L;
	public static final double[] EQUAL_FREQ = new double[]{0.25, 0.25, 0.25, 0.25};
	public static final int DIMENSION = 4;
	public static final Bounds<Double> SPECTRA_BOUNDS = new DefaultBounds(1.0, 0.0, DIMENSION);
	
	
	public SpectraParameter2(SpectrumType type){
		this(type.getCode());
	}
	
	public SpectraParameter2(int type){
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
				freq[i] = 1+MathUtils.nextInt(100)+MathUtils.nextInt(100);
				sum += freq[i];
			}
			for (int i = 0; i < freq.length; i++) {
				freq[i] /= sum;
				if(freq[i]<=0){
					System.err.println(Arrays.toString(freq) +"\t"+ sum);
				}
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

    public SpectraParameter2(double[] frequencies) {
//    	super(frequencies);
    	this.values = new double[values.length];
        System.arraycopy(values, 0, this.values, 0, values.length);
    	setId("spectra");
    	

        double sum = getSumOfFrequencies(frequencies);
    	if(values.length!=DIMENSION){
    		throw new IllegalArgumentException("Frequencies should have 4 elements, frequencies.length= "+getDimension());
    	}
        if (Math.abs(sum - 1.0) > 1e-8) {
            throw new IllegalArgumentException("Frequencies do not sum to 1, they sum to " + sum);
        }
    	
		addBounds(SPECTRA_BOUNDS);
		if(!isWithinBounds()){
			throw new IllegalArgumentException("Frequencies out of bounds 0 < f < 1\t"+ Arrays.toString(frequencies)); 
		}

    }

    private static double getSumOfFrequencies(double[] frequencies) {
        double total = 0.0;
        for (int i = 0; i < DIMENSION; i++) {
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
        for (int i = 1; i < DIMENSION; i++) {
            frequencies[i] += frequencies[i - 1];
        }
        return frequencies;
    }
    

    public String diagnostic(){
    	String diag = Arrays.toString(getFrequencies());
    	return diag;
    }
    
//    private void addBounds(){
//		addBounds(SPECTRA_BOUNDS);
//	}

	public double getStoredFrequency(int state) {
//		storedValues
		
		return storedFrequency[state];
	}
	private double[] storedFrequency = new double[DIMENSION];

}
///////////////
	abstract class AbstractSpectra extends Parameter.Abstract{
		
	
    
	private static final long serialVersionUID = 8795190465907061489L;

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

    public final int getDimension() {
        return values.length;
    }

    public final int getSize() {
        return getDimension();
    }

    public final double getParameterValue(int i) {
        return values[i];
    }

    /**
     * Defensively returns copy of parameter array.
     *
     * @return a copy of the parameter values
     */
    public final double[] getParameterValues() {

        double[] copyOfValues = new double[values.length];
        System.arraycopy(values, 0, copyOfValues, 0, copyOfValues.length);
        return copyOfValues;
    }

    /**
     * Do not write to the returned array directly!!
     *
     * @return the parameter values
     */
//    public final double[] inspectParameterValues() {
//        return values;
//    }

    public Bounds<Double> getBounds() {
        if (bounds == null) {
            throw new NullPointerException(getParameterName() + " parameter: Bounds not set");
        }
        return bounds;
    }

    public String getParameterName() {
        return getId();
    }

    //********************************************************************
    // SETTERS
    //********************************************************************

    /**
     * Can only be called before store is called. If it results in new
     * dimensions, then the value of the first dimension is copied into the new dimensions.
     */
    public void setDimension(int dim) {
    	throw new IllegalArgumentException("Can not setDimension");
//        final int oldDim = getDimension();
//        if (oldDim == dim) {
//            return;
//        }
//
//        assert storedValues == null :
//                "Can't change dimension after store has been called! storedValues=" +
//                        Arrays.toString(storedValues) + " bounds=" + bounds;
//
//
//        double[] newValues = new double[dim];
//        // copy over new values, min in case new dim is smaller
//        System.arraycopy(values, 0, newValues, 0, Math.min(oldDim, dim));
//        // fill new values with first item
//        for (int i = oldDim; i < dim; i++) {
//            newValues[i] = values[0];
//        }
//        values = newValues;
//
//        if (bounds != null) {
//            //assert oldDim < dim :  "Can't decrease dimension when bounds are set";
//            for (int k = 1; k < oldDim; ++k) {
//                assert ((double) bounds.getLowerLimit(k) == bounds.getLowerLimit(0)) &&
//                        ((double) bounds.getUpperLimit(k) == bounds.getUpperLimit(0)) :
//                        "Can't change dimension when bounds are not all equal";
//            }
//            final double low = bounds.getLowerLimit(0);
//            final double high = bounds.getUpperLimit(0);
//            bounds = null;
//            addBounds(low, high);
//        }
    }

    public void addDimension(int index, double value) {
    	throw new IllegalArgumentException("Can not addDimension");
    }
    public double removeDimension(int index) {
    	throw new IllegalArgumentException("Can not removeDimension");
    }


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
    public void setParameterValueQuietly(int dim, double value) {
        values[dim] = value;
    }


    /**
     * Sets the values of the parameter and notify that all values of the parameter have changed.
     *
     * @param i   index of the value
     * @param val to value to set
     */
    public void setParameterValueNotifyChangedAll(int i, double val) {
        values[i] = val;
        fireParameterChangedEvent(i, Parameter.ChangeType.ALL_VALUES_CHANGED);
    }

    protected final void storeValues() {
        // no need to pay a price in a very common call for one-time rare usage
        //hasBeenStored = true;
        if (storedValues == null || storedValues.length != values.length) {
            storedValues = new double[values.length];
        }
        System.arraycopy(values, 0, storedValues, 0, storedValues.length);
    }

    protected final void restoreValues() {

        //swap the arrays
        double[] temp = storedValues;
        storedValues = values;
        values = temp;

        //if (storedValues != null) {
        //	System.arraycopy(storedValues, 0, values, 0, values.length);
        //} else throw new RuntimeException("restore called before store!");
    }

    /**
     * Nothing to do
     */
    protected final void acceptValues() {
    }

    protected final void adoptValues(Parameter source) {
        // todo bug ? bounds not adopted?

        if (getDimension() != source.getDimension()) {
            throw new RuntimeException("The two parameters don't have the same number of dimensions");
        }

        for (int i = 0, n = getDimension(); i < n; i++) {
            values[i] = source.getParameterValue(i);
        }
    }

    protected double[] values;

    protected double[] storedValues;

    // same as !storedValues && !bounds
    //private boolean hasBeenStored = false;
    private Bounds<Double> bounds = null;

//    public void addBounds(double lower, double upper) {
//        addBounds(new DefaultBounds(upper, lower, getDimension()));
//    }
}