package srp.spectrum.operator;

import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.util.FastMath;

import srp.spectrum.SpectraParameter;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;

import com.google.common.primitives.Ints;

import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class SwapSubColumnSpectrumOperator extends AbstractSwapSpectrumOperator {

	public static final String OPERATOR_NAME = SwapSubColumnSpectrumOperator.class.getSimpleName();
	public static final SpectrumOperation OP = SpectrumOperation.SWAP_SUBCOLUMN;
//    private Parameter parameter = null;
    private final int[] parameterWeights;
//    private double swapHapCount = 0.05;


	private int swapSpectrumCount;
	private double autoOptimize;
//	private double scaleFactor;

	
	public SwapSubColumnSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			int swapSpectrumCount, CoercionMode mode) {
		super(spectrumModel, mode, true);
		
		this.swapSpectrumCount = swapSpectrumCount;
//		this.spectrumCount = spectrumModel.getSpectrumCount();
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }
		convertToAutoOptimize(this.swapSpectrumCount);


//		scaleFactor = (int) (spectrumLength*0.01);
//		if (scaleFactor <1) {
//			scaleFactor = 1;
//		}
		

	}	
	
	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();
		int siteIndex = MathUtils.nextInt(spectrumLength);

		int[] spectrumIndexs = randomSiteHashSet(swapSpectrumCount, spectrumLength);

		for (int i = 0; i < swapSpectrumCount; i++) {
			int s = spectrumIndexs[i];
//			SpectraParameter spectra = spectrum.getSpectra(spectrumIndexs[i]);
			SpectraParameter spectra = spectrumModel.getSpectrum(s).getSpectra(siteIndex);
			swapFrequency(spectra);
		}

		spectrumModel.setSpectrumOperationRecord(OP, spectrumIndexs, siteIndex);
		spectrumModel.endSpectrumOperation();

		return 0.0;
	}



	@Override
	public void setCoercableParameter(double autoOpt) {
		convertFromAutoOptimizeToValue(autoOpt);
	}
	@Override
	public double getCoercableParameter() {
	    return autoOptimize;
	}


	private void convertFromAutoOptimizeToValue(double autoOpt) {
    	autoOptimize = autoOpt;
//    	swapHapCount =  1 + (int) Math.exp(autoOptimize*scaleFactor);
    	swapSpectrumCount =  1 + (int) FastMath.exp(autoOptimize);
		
		checkParameterIsValid();
    }

	private double convertToAutoOptimize(int length) {
		swapSpectrumCount = length;
		checkParameterIsValid();
		autoOptimize = Math.log(swapSpectrumCount - 1);
//		autoOptimize = Math.log(swapHapCount - 1)/scaleFactor;
	    return autoOptimize;
	}

	private void checkParameterIsValid() {
		if (swapSpectrumCount > spectrumLength){
			swapSpectrumCount = spectrumLength;
		}
	}
	


    @Override
	public double getRawParameter() {
        return swapSpectrumCount;
    }


    @Override
	public final String getPerformanceSuggestion() {

    	String s = "Tuning "+swapSpectrumCount; 
    	return s;


    }

    @Override
	public String toString() {
        return getOperatorName() + "(windowsize=" + swapSpectrumCount + ")";
    }


	@Override
	public SpectrumOperation getSpectrumOperation() {
		return OP;
	}


	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
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
