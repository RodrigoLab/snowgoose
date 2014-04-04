package srp.evolution.spectrum;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;

public abstract class AbstractSpectra extends Parameter.Abstract {

	
    /**
	 * 
	 */
	private static final long serialVersionUID = -1845301617682337277L;
	
//	protected double[] values;
//    protected double[] storedValues;
    
    protected int dimension;
	Bounds<Double> bounds = null;
    
	
//	public AbstractSpectra() {
////		super();
//	}
//
//	public AbstractSpectra(String name) {
//		super(name);
//	}
	
	public abstract double[] getFrequencies();

	public abstract double getFrequency(int f);

	protected abstract void storeValues();
	protected abstract void restoreValues();
	protected abstract void acceptValues();
	

	@Override
	public int getDimension() {
	    return dimension;
	}

	@Override
	public int getSize() {
	    return dimension;
	}
	
	@Override
	public void addBounds(Bounds<Double> boundary) {
	    if (bounds == null) {
	        bounds = boundary;
	    } else {
	        throw new IllegalArgumentException("Should not call addBounds twice");
	        // can't change dimension after bounds are added!
	    }
	
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

	//IllegalArgumentException below!!
	@Override
	public final double getParameterValue(int i) {
		throw new IllegalArgumentException("Use getCategory/Frequency");
	}

	@Override
	public final double[] getParameterValues() {
		throw new IllegalArgumentException("Use getCategory/Frequency");
	}

	@Override
	public final void setParameterValue(int cat, double val) {
	    	throw new IllegalArgumentException("Use setCategory/Frequency");
	//        fireParameterChangedEvent(cat, Parameter.ChangeType.VALUE_CHANGED);
	    }

	/**
	 * Sets the value of the parameter without firing a changed event.
	 *
	 * @param cat   the index of the parameter dimension
	 * @param value the value to set
	 */
	@Override
	public final void setParameterValueQuietly(int cat, double value) {
		throw new IllegalArgumentException("Use setCategory/Frequency");
	}

	/**
	 * Sets the values of the parameter and notify that all values of the parameter have changed.
	 *
	 * @param i   index of the value
	 * @param val to value to set
	 */
	@Override
	public final void setParameterValueNotifyChangedAll(int i, double val) {
	    	throw new IllegalArgumentException("Use getCategory/Frequency");
	//    	fireParameterChangedEvent(i, Parameter.ChangeType.ALL_VALUES_CHANGED);
	    }

	@Override
	public final void setDimension(int dim) {
		throw new IllegalArgumentException("Can not setDimension");
	}

	@Override
	public final void addDimension(int index, double value) {
		throw new IllegalArgumentException("Can not addDimension");
	}

	@Override
	public final double removeDimension(int index) {
		throw new IllegalArgumentException("Can not removeDimension");
	}

	@Override
	public final void adoptValues(Parameter source) {
		throw new IllegalArgumentException("Can not adoptValues");
	}

}