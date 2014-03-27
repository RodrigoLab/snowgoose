package srp.spectrum;

import java.util.Arrays;

import srp.dr.evolution.datatype.ShortReads;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

public class CategorySpectraParameter extends AbstractSpectra{
//public class CategorySpectraParameter extends SpectraParameter{
	
	private static final long serialVersionUID = 3708765314559863330L;

	public static final int TOTAL_STATE_COUNT = ShortReads.INSTANCE.getAmbiguousStateCount();

	public static final int ONE_CAT = 4;
	public static final int TWO_CATS = 10;
	public static final int ALL_CATS = 10;

	public static int DIMENSION = 1;
	public static final Bounds<Double> CATEGORY_BOUNDS = new DefaultBounds(9.0, 0.0, DIMENSION);
	
    protected int[] values;
    protected int[] storedValues;
			
	public CategorySpectraParameter() {
		setId("spectra_category");
		
		addBounds(CATEGORY_BOUNDS);
		this.dimension = DIMENSION;
		this.values = new int[DIMENSION];
	    this.storedValues = new int[DIMENSION];
	}

	public CategorySpectraParameter(CategoryType type) {
		this();
		int cat;
		switch (type) {
		case SINGLE:
			cat = MathUtils.nextInt(ONE_CAT);
			break;
		case TWOWAYS:
			cat = MathUtils.nextInt(TWO_CATS);
			break;
		case THREEWAYS:		
		case RANDOM:
			cat = MathUtils.nextInt(ALL_CATS);
			break;
		default:
			throw new IllegalArgumentException("Invalid type: "+type);
		}
		setCategory(cat);
	}
	
    public CategorySpectraParameter(int cat) {
    	this();
    	setCategory(cat);
		if(!isWithinBounds()){
			throw new IllegalArgumentException("Category out of bounds 0 < f < "+ CATEGORY_BOUNDS.getUpperLimit(0)+"\t"+ cat); 
		}

    }
    public void setCategory(int cat) {
//    	setParameterValue(i, 0);
    	values[0] = cat;
    	fireParameterChangedEvent(cat, Parameter.ChangeType.ALL_VALUES_CHANGED);
    }
    
    public int getCategory() {
    	return values[0];
    }
    
    public double[] getFrequencies() {
    	return CATEGORIES[values[0]];
    }
    
    public double getFrequency(int f) {
    	return CATEGORIES[values[0]][f];
    }

	
	@Override
	public boolean isWithinBounds() {
	    Bounds<Double> bounds = getBounds();
	    final double value = getCategory();
	    if (value < bounds.getLowerLimit(0) || value > bounds.getUpperLimit(0)) {
	    	return false;
	    }
	    return true;
	}


	public enum CategoryType{
		SINGLE,
		TWOWAYS,
		THREEWAYS,		
		RANDOM, 
	}
	
	public static final double[][] CATEGORIES= new double[][]{
		{0.97, 0.01, 0.01, 0.01},
		{0.01, 0.97, 0.01, 0.01},
		{0.01, 0.01, 0.97, 0.01},
		{0.01, 0.01, 0.01, 0.97},
		
		{0.49, 0.49, 0.01, 0.01},
		{0.49, 0.01, 0.49, 0.01},
		{0.49, 0.01, 0.01, 0.49},
		{0.01, 0.49, 0.49, 0.01},
		{0.01, 0.49, 0.01, 0.49},
		{0.01, 0.01, 0.49, 0.49},
	};

    @Override
	protected final void storeValues() {
//        System.arraycopy(values, 0, storedValues, 0, DIMENSION);
        storedValues[0] = values[0];
    }

    @Override
	protected final void restoreValues() {
        //swap the arrays
//        int[] temp = storedValues;
//        storedValues = values;
//        values = temp;
    	values[0] = storedValues[0]; 
    }

    /**
     * Nothing to do
     */
    @Override
	protected final void acceptValues() {
    }



}
