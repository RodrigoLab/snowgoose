package srp.spectrum.operator;

import srp.spectrum.SpectraParameter;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class SwapColumnSpectrumOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = SwapColumnSpectrumOperator.class.getSimpleName();
	public static final SpectrumOperation OP = SpectrumOperation.SWAP_COLUMN;
//    private Parameter parameter = null;
    private final int[] parameterWeights;
    private double delta = 0.05;
    
	public SwapColumnSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			double delta, CoercionMode mode) {
		super(spectrumModel, mode);
		
		
		this.delta = delta;
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }

	}

	private double[] debugList = new double[8];
	
	
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();
		int siteIndex = MathUtils.nextInt(spectrumLength);

//        SpectraParameter[] spectra = new SpectraParameter[spectrumCount];
//        double[] d = new double[spectrumCount];
//        double[] scalar1 = new double[spectrumCount];
//        double[] scalar2 = new double[spectrumCount];
//        int[] dim1 = new int[spectrumCount];
//		int[] dim2 = new int[spectrumCount];
		
        
		
		for (int i = 0; i < spectrumCount; i++) {
			SpectraParameter spectra = spectrumModel.getSpectrum(i).getSpectra(siteIndex);
			
			// get two dimensions
			
//	        int dim1 = MathUtils.nextInt(DIMENSION);
//	        int dim2;// = dim1;
//	        do {
//	            dim2 = MathUtils.nextInt(DIMENSION);
//	        }while (dim1 == dim2);
//	        
//	        double scalar1 = parameter.getParameterValue(dim1);
//	        double scalar2 = parameter.getParameterValue(dim2);
//			
			int dim1 = 0;
	        double scalar1=0;
	        for (int j = 0; j < DIMENSION; j++) {
	        	scalar1 = spectra.getParameterValue(j);
	        	if(scalar1==1){
	        		dim1=j;
	        		break;
	        	}
			}
	        int dim2;// = dim1;
	        do {
	            dim2 = MathUtils.nextInt(DIMENSION);
	        }while (dim1 == dim2);
	        double scalar2 = spectra.getParameterValue(dim2);


	        spectra.setParameterValue(dim1, scalar2);
	        spectra.setParameterValue(dim2, scalar1);

//			
//	        // get two dimensions
//	        dim1[i] = MathUtils.nextInt(DIMENSION);
//	        do {
//	            dim2[i] = MathUtils.nextInt(DIMENSION);
//	        }while (dim1[i] == dim2[i]);
//    
//	        scalar1[i] = spectra[i].getParameterValue(dim1[i]);
//	        scalar2[i]= spectra[i].getParameterValue(dim2[i]);
//	
//	        d[i] = MathUtils.nextDouble() * delta;
//	        d[i] = delta;
//	        scalar1[i] -= d[i];
//	        scalar2[i] += d[i];
//	
//	        
//	
//	        if (scalar1[i] < BOUNDS_LOWER ||
//	                scalar1[i] > BOUNDS_UPPER ||
//	                scalar2[i] < BOUNDS_LOWER ||
//	                scalar2[i] > BOUNDS_UPPER ) {
////	        	System.err.println("throw");
//	            throw new OperatorFailedException("proposed values out of range!");
//	        }
	        
		}
//		for (int i = 0; i < spectrumCount; i++) {
//			
//			spectra[i].setParameterValue(dim1[i], scalar1[i]);
//			spectra[i].setParameterValue(dim2[i], scalar2[i]);
//
//		}
        // symmetrical move so return a zero hasting ratio
		spectrumModel.setSpectrumOperationRecord(OP, siteIndex);
		
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
//		String s = record.getSpectrumIndex() +"\t"+ record.getColumnIndex() +"\n";
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
