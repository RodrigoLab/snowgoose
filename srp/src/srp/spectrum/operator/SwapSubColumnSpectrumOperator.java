package srp.spectrum.operator;

import java.util.HashSet;
import java.util.Set;

import srp.spectrum.SpectraParameter;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;

import com.google.common.primitives.Ints;

import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class SwapSubColumnSpectrumOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = SwapSubColumnSpectrumOperator.class.getSimpleName();
	public static final SpectrumOperation OP = SpectrumOperation.SWAP_SUBCOLUMN;
//    private Parameter parameter = null;
    private final int[] parameterWeights;
//    private double swapHapCount = 0.05;


	private double[] debugList = new double[8];
	private int swapSpectrumCount;
	private double autoOptimize;
	private int spectrumCount;
//	private double scaleFactor;

	
	public SwapSubColumnSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			int swapHapCount, CoercionMode mode) {
		super(spectrumModel, mode);
		
		
		this.swapSpectrumCount = swapHapCount;
		this.spectrumCount = spectrumModel.getSpectrumCount();
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

//        SpectraParameter[] spectra = new SpectraParameter[spectrumCount];
//        double[] d = new double[spectrumCount];
//        double[] scalar1 = new double[spectrumCount];
//        double[] scalar2 = new double[spectrumCount];
//        int[] dim1 = new int[spectrumCount];
//		int[] dim2 = new int[spectrumCount];
		
		
		Set<Integer> generated = new HashSet<Integer>();
		while (generated.size() < swapSpectrumCount)
		{
		    Integer next = MathUtils.nextInt(spectrumCount);
		    generated.add(next);
		}
//		generated.toArray(siteIndex);
		int[] spectrumIndexs = Ints.toArray(generated);
		
		//		System.out.println(Arrays.toString(siteIndex));
		for (int h = 0; h < swapSpectrumCount; h++) {
			int i = spectrumIndexs[h];
			SpectraParameter spectra = spectrumModel.getSpectrum(i).getSpectra(siteIndex);

//			// get any two dims and swap
//	        int dim1 = MathUtils.nextInt(DIMENSION);
//	        int dim2;// = dim1;
//	        do {
//	            dim2 = MathUtils.nextInt(DIMENSION);
//	        }while (dim1 == dim2);
//	        
//	        double scalar1 = parameter.getParameterValue(dim1);
//	        double scalar2 = parameter.getParameterValue(dim2);
	        
	        
			// get freq==1 and swap with others
			int dim1 = 0;
			double scalar1 = 0;
			for (int d = 0; d < DIMENSION; d++) {
				scalar1 = spectra.getParameterValue(d);
				if (scalar1 == 1) {
					dim1 = d;
					break;
				}
			}
			int dim2;// = dim1;
			do {
				dim2 = MathUtils.nextInt(DIMENSION);
			} while (dim1 == dim2);
			double scalar2 = spectra.getParameterValue(dim2);        

	        spectra.setParameterValue(dim1, scalar2);
	        spectra.setParameterValue(dim2, scalar1);


	        
		}

		spectrumModel.setSpectrumOperationRecord(OP, spectrumIndexs, siteIndex);
		//TODO: finish implement calculation/store/restore
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
    	swapSpectrumCount =  1 + (int) Math.exp(autoOptimize);
		
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
//		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
//		String s = record.getSpectrumIndex() +"\t"+ record.getColumnIndex() +"\n";
//		s+= spectrumModel.diagnostic() +"\n";
//		s += Arrays.toString(debugList);
//		spectrumModel.restoreModelState();
    	String s = "Tuning "+swapSpectrumCount; 
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
