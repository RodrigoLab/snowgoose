package srp.spectrum.operator;

import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class SwapSingleSpectrumOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = SwapSingleSpectrumOperator.class.getSimpleName();
	public static final SpectrumOperation OP = SpectrumOperation.DELTA_SINGLE;
//    private Parameter parameter = null;
    private final int[] parameterWeights;
    private int[] siteIndex;
    private double delta = Double.NaN;

    
	public SwapSingleSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			CoercionMode mode) {
		super(spectrumModel, mode);
		
		
//		this.delta = delta;
		this.siteIndex = new int[1];
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }

	}

	private double[] debugList = new double[8];
	private int mm;
	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();

		int spectrumIndex = MathUtils.nextInt(spectrumCount);
		siteIndex[0] = MathUtils.nextInt(spectrumLength);
		
		Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
		SpectraParameter parameter = spectrum.getSpectra(siteIndex[0]);
        // get two dimensions
		
        int dim1 = MathUtils.nextInt(DIMENSION);
        int dim2;// = dim1;
        do {
            dim2 = MathUtils.nextInt(DIMENSION);
        }while (dim1 == dim2);
        
        double scalar1 = parameter.getParameterValue(dim1);
        double scalar2 = parameter.getParameterValue(dim2);
        
//		int dim1 = 0;
//        double scalar1=0;
//        for (int i = 0; i < DIMENSION; i++) {
//        	scalar1 = parameter.getParameterValue(i);
//        	if(scalar1==1){
//        		dim1=i;
//        		break;
//        	}
//		}
//
//        int dim2;// = dim1;
//        do {
//            dim2 = MathUtils.nextInt(DIMENSION);
//        }while (dim1 == dim2);
//        double scalar2 = parameter.getParameterValue(dim2);
        
        
//        double d = MathUtils.nextDouble() * delta;
//        d = delta;
//        scalar1 -= d;
//        scalar2 += d;



        parameter.setParameterValue(dim1, scalar2);
        parameter.setParameterValue(dim2, scalar1);

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






