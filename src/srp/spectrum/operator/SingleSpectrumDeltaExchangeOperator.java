package srp.spectrum.operator;

import java.util.Arrays;

import srp.haplotypes.HaplotypeModel;
import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.OperatorUtils;
import dr.math.MathUtils;

public class SingleSpectrumDeltaExchangeOperator extends AbstractSpectrumOperator {

	public final static String OPERATOR_NAME = SingleSpectrumDeltaExchangeOperator.class.getSimpleName();

//    private Parameter parameter = null;
    private final int[] parameterWeights;
    private double delta = 0.02;
    private boolean isIntegerOperator = false;
    
	public SingleSpectrumDeltaExchangeOperator(SpectrumAlignmentModel spectrumModel, 
			double delta, CoercionMode mode) {
		super(spectrumModel, mode);
		
		
		this.delta = delta;
        setWeight(1.0);
        this.isIntegerOperator = false;

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }

	}

	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();

//		int[] posChar = alignmentMapping.getNextBase();
//		spectrumModel.swapHaplotypeSingleBase(OP, posChar);
		int spectrumIndex = MathUtils.nextInt(spectrumCount);
		int site = MathUtils.nextInt(spectrumLength);
		
		System.err.println(spectrumIndex +"\t"+ site);
		
		Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
		SpectraParameter parameter = spectrum.getSpecturm(site);
        // get two dimensions
        final int dim = parameter.getDimension();
        final int dim1 = MathUtils.nextInt(dim);
        int dim2;// = dim1;
        do {
            dim2 = MathUtils.nextInt(dim);
        }while (dim1 == dim2);

        double scalar1 = parameter.getParameterValue(dim1);
        double scalar2 = parameter.getParameterValue(dim2);

//        if (isIntegerOperator) {
//            int d = MathUtils.nextInt((int) Math.round(delta)) + 1;
//
//            if (parameterWeights[dim1] != parameterWeights[dim2]) throw new RuntimeException();
//            scalar1 = Math.round(scalar1 - d);
//            scalar2 = Math.round(scalar2 + d);
//        } else {

            // exchange a random delta
            final double d = MathUtils.nextDouble() * delta;
            scalar1 -= d;
//            if (parameterWeights[dim1] != parameterWeights[dim2]) {
//                scalar2 += d * (double) parameterWeights[dim1] / (double) parameterWeights[dim2];
//            } else {
                scalar2 += d;
//            }

//        }
        Bounds<Double> bounds = parameter.getBounds();

        if (scalar1 < bounds.getLowerLimit(dim1) ||
                scalar1 > bounds.getUpperLimit(dim1) ||
                scalar2 < bounds.getLowerLimit(dim2) ||
                scalar2 > bounds.getUpperLimit(dim2)) {
            throw new OperatorFailedException("proposed values out of range!");
        }
        System.out.println(scalar1 +"\t"+ scalar2);
        parameter.setParameterValue(dim1, scalar1);
        parameter.setParameterValue(dim2, scalar2);

        // symmetrical move so return a zero hasting ratio
//        spectrum = spectrumModel.getSpectrum(spectrumIndex);
//        parameter = spectrum.getSpecturm(site);
//		System.out.println(Arrays.toString(parameter.getFrequencies()));
		
		
		
		
		
		spectrumModel.endSpectrumOperation();

		return 0.0;
	}

	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}


    public double getCoercableParameter() {
        return Math.log(delta);
    }

    public void setCoercableParameter(double value) {
        delta = Math.exp(value);
    }

    public double getRawParameter() {
        return delta;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public final String getPerformanceSuggestion() {

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();

//        double d = OperatorUtils.optimizeWindowSize(delta, parameter.getParameterValue(0) * 2.0, prob, targetProb);
        double d = OperatorUtils.optimizeWindowSize(delta, 0.25 , prob, targetProb);

        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try decreasing delta to about " + d;
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try increasing delta to about " + d;
        } else return "";
    }

    public String toString() {
        return getOperatorName() + "(windowsize=" + delta + ")";
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
