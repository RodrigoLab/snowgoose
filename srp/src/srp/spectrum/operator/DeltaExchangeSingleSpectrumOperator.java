package srp.spectrum.operator;

import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class DeltaExchangeSingleSpectrumOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = DeltaExchangeSingleSpectrumOperator.class.getSimpleName();
	public static final SpectrumOperation OP = SpectrumOperation.DELTA_SINGLE;
//    private Parameter parameter = null;
    private final int[] parameterWeights;
    
    private int[] siteIndex;
    private double delta = 0.1;

    private double[] debugList = new double[8];
	private int mm;
//	public int failcount = 0;
	
	public DeltaExchangeSingleSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			double delta, CoercionMode mode) {
		super(spectrumModel, mode);
		
		
		this.delta = delta;
		this.siteIndex = new int[1];
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }

	}
	
	
	@Override
	public double doOperation() throws OperatorFailedException {
	
		spectrumModel.startSpectrumOperation();

		int spectrumIndex = MathUtils.nextInt(spectrumCount);
		siteIndex[0] = MathUtils.nextInt(spectrumLength);
		
		Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
		SpectraParameter parameter = spectrum.getSpectra(siteIndex[0]);
        // get two dimensions

        final int dim1 = MathUtils.nextInt(DIMENSION);
        int dim2;// = dim1;
        do {
            dim2 = MathUtils.nextInt(DIMENSION);
        }while (dim1 == dim2);

        double scalar1 = parameter.getParameterValue(dim1);
        double scalar2 = parameter.getParameterValue(dim2);

//		double d = delta;
        double d = MathUtils.nextDouble() * delta;
//        d = delta;
        scalar1 -= d;
        scalar2 += d;


    	if (scalar1 < BOUNDS_LOWER ||
                scalar1 > BOUNDS_UPPER ||
                scalar2 < BOUNDS_LOWER ||
                scalar2 > BOUNDS_UPPER ) {        	
//    		failcount ++;
//    		System.out.println(d +"\t"+ spectrumIndex +"\t"+ siteIndex[0] +"\t"+ scalar1 +"\t"+ scalar2 +"\t"+ Arrays.toString(parameter.getFrequencies()));
            throw new OperatorFailedException("proposed values out of range!");
        }

        parameter.setParameterValue(dim1, scalar1);
        parameter.setParameterValue(dim2, scalar2);

        // symmetrical move so return a zero hasting ratio
		spectrumModel.setSpectrumOperationRecord(OP, spectrumIndex, siteIndex);
		
		spectrumModel.endSpectrumOperation();

		return 0.0;
	}

	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}


    @Override
	public double getCoercableParameter() {
//    	double t = Math.log(delta/(1-delta));
//    	return t;
//        return Math.log(1.0 / delta - 1.0);
        return Math.log(delta);
    }

    @Override
	public void setCoercableParameter(double value) {
//    	mm++;
//    	if(mm%1000 == 0){
//    		System.out.println(value +"\t"+ delta +"\t"+ getAcceptanceProbability());
//    	}
        delta = Math.exp(value);
//        double t = Math.exp(value);
//        delta = t/(t+1);

    }

    @Override
	public double getRawParameter() {
        return delta;
    }

    @Override
	public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    @Override
	public final String getPerformanceSuggestion() {
//		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
//		String s = record.getSpectrumIndex() +"\t"+ Arrays.toString(record.getAllSiteIndexs()) +"\n";
//		s+= spectrumModel.diagnostic() +"\n";
//		s += Arrays.toString(debugList);
//		spectrumModel.restoreModelState();
    	String s = "Tuning "+delta; 
    	return s;
    	
//        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
//        double targetProb = getTargetAcceptanceProbability();
//
////        double d = OperatorUtils.optimizeWindowSize(delta, parameter.getParameterValue(0) * 2.0, prob, targetProb);
//        double d = OperatorUtils.optimizeWindowSize(delta, 0.25 , prob, targetProb);
//
//        if (prob < getMinimumGoodAcceptanceLevel()) {
//            return "Try decreasing delta to about " + d;
//        } else if (prob > getMaximumGoodAcceptanceLevel()) {
//            return "Try increasing delta to about " + d;
//        } else return "";
    }

    @Override
	public String toString() {
        return getOperatorName() + "(windowsize=" + delta + ")";
    }

	@Override
	public SpectrumOperation getSpectrumOperation() {
		return OP;
	}

	
}







////////////////////////

class DeltaExchangeOperator {//extends AbstractCoercableOperator {

    public DeltaExchangeOperator(Parameter parameter, double delta) {

//        super(CoercionMode.COERCION_ON);

        this.parameter = parameter;
        this.delta = delta;
//        setWeight(1.0);
        this.isIntegerOperator = false;

        parameterWeights = new int[parameter.getDimension()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }
    }

    public DeltaExchangeOperator(Parameter parameter, int[] parameterWeights, double delta, double weight, boolean isIntegerOperator, CoercionMode mode) {

//        super(mode);

        this.parameter = parameter;
        this.delta = delta;
//        setWeight(weight);
        this.isIntegerOperator = isIntegerOperator;
        this.parameterWeights = parameterWeights;

        if (isIntegerOperator && delta != Math.round(delta)) {
            throw new IllegalArgumentException("Can't be an integer operator if delta is not integer");
        }
    }

  
    /**
     * change the parameter and return the hastings ratio.
     * performs a delta exchange operation between two scalars in the vector
     * and return the hastings ratio.
     */
    public final double doOperation() throws OperatorFailedException {

        // get two dimensions
        final int dim = parameter.getDimension();
        final int dim1 = MathUtils.nextInt(dim);
        int dim2 = dim1;
        while (dim1 == dim2) {
            dim2 = MathUtils.nextInt(dim);
        }

        double scalar1 = parameter.getParameterValue(dim1);
        double scalar2 = parameter.getParameterValue(dim2);

        if (isIntegerOperator) {
            int d = MathUtils.nextInt((int) Math.round(delta)) + 1;

            if (parameterWeights[dim1] != parameterWeights[dim2]) throw new RuntimeException();
            scalar1 = Math.round(scalar1 - d);
            scalar2 = Math.round(scalar2 + d);
        } else {

            // exchange a random delta
            final double d = MathUtils.nextDouble() * delta;
            scalar1 -= d;
            if (parameterWeights[dim1] != parameterWeights[dim2]) {
                scalar2 += d * (double) parameterWeights[dim1] / (double) parameterWeights[dim2];
            } else {
                scalar2 += d;
            }

        }
        Bounds<Double> bounds = parameter.getBounds();

        if (scalar1 < bounds.getLowerLimit(dim1) ||
                scalar1 > bounds.getUpperLimit(dim1) ||
                scalar2 < bounds.getLowerLimit(dim2) ||
                scalar2 > bounds.getUpperLimit(dim2)) {
            throw new OperatorFailedException("proposed values out of range!");
        }
        parameter.setParameterValue(dim1, scalar1);
        parameter.setParameterValue(dim2, scalar2);

        // symmetrical move so return a zero hasting ratio
        return 0.0;
    }
//
//    // Interface MCMCOperator
//    public final String getOperatorName() {
//        return parameter.getParameterName();
//    }
//
//    public double getCoercableParameter() {
//        return Math.log(delta);
//    }
//
//    public void setCoercableParameter(double value) {
//        delta = Math.exp(value);
//    }
//
//    public double getRawParameter() {
//        return delta;
//    }
//
//    public double getTargetAcceptanceProbability() {
//        return 0.234;
//    }
//
//    public final String getPerformanceSuggestion() {
//
//        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
//        double targetProb = getTargetAcceptanceProbability();
//
//        double d = OperatorUtils.optimizeWindowSize(delta, parameter.getParameterValue(0) * 2.0, prob, targetProb);
//
//
//        if (prob < getMinimumGoodAcceptanceLevel()) {
//            return "Try decreasing delta to about " + d;
//        } else if (prob > getMaximumGoodAcceptanceLevel()) {
//            return "Try increasing delta to about " + d;
//        } else return "";
//    }
//
//    public String toString() {
//        return getOperatorName() + "(windowsize=" + delta + ")";
//    }

    // Private instance variables

    
    private Parameter parameter = null;
    private final int[] parameterWeights;
    private double delta = 0.02;
    private boolean isIntegerOperator = false;
}
