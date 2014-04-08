package srp.operator.spectrum;

import org.apache.commons.math3.util.FastMath;

import srp.evolution.OperationType;
import srp.evolution.spectrum.SpectraParameter;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class SwapColumnSpectrumOperator extends AbstractSwapSpectrumOperator {

	public static final String OPERATOR_NAME = SwapColumnSpectrumOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.COLUMN;
	
    private final int[] parameterWeights;
    private double delta;
    
    public SwapColumnSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			CoercionMode mode) {
		this(spectrumModel, mode, true);
    }
		
	public SwapColumnSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			CoercionMode mode, boolean random) {
		super(spectrumModel, mode, random);

		setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }
        

	}

	
	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startAlignmentModelOperation();
		int siteIndex = MathUtils.nextInt(spectrumLength);
		
		for (int i = 0; i < spectrumCount; i++) {
			SpectraParameter spectra = spectrumModel.getSpectrum(i).getSpectra(siteIndex);
			swapFrequency(spectra);
		}

        // symmetrical move so return a zero hasting ratio
		spectrumModel.setOperationRecord(OP, fixSpectrumIndexArray, siteIndex);
		
		spectrumModel.endAlignmentModelOperation();

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
        delta = FastMath.exp(value);
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
	public OperationType getSpectrumOperation() {
		return OP;
	}

	
}


//
//public double getMinimumAcceptanceLevel() {
//    return 0.05;
//}
//
//public double getMaximumAcceptanceLevel() {
//    return 0.50;
//}
//
//public double getMinimumGoodAcceptanceLevel() {
//    return 0.10;
//}
//
//public double getMaximumGoodAcceptanceLevel() {
//    return 0.40;
//}
