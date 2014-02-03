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

    
	public SwapSingleSpectrumOperator(SpectrumAlignmentModel spectrumModel) {
		super(spectrumModel, CoercionMode.COERCION_OFF);
		
		
//		this.delta = delta;
		this.siteIndex = new int[1];
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }

	}

	private double[] debugList = new double[8];

	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();

		int spectrumIndex = MathUtils.nextInt(spectrumCount);
		siteIndex[0] = MathUtils.nextInt(spectrumLength);
		
		Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
		SpectraParameter spectra = spectrum.getSpectra(siteIndex[0]);

//		// get any two dims and swap
//        int dim1 = MathUtils.nextInt(DIMENSION);
//        int dim2;// = dim1;
//        do {
//            dim2 = MathUtils.nextInt(DIMENSION);
//        }while (dim1 == dim2);
//        
//        double scalar1 = parameter.getParameterValue(dim1);
//        double scalar2 = parameter.getParameterValue(dim2);
        
        
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
        return Double.NaN;
    }

    @Override
	public void setCoercableParameter(double value) {
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
    	String s = "No tuning"; 
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
        return getOperatorName() + "(Single)";
    }

	@Override
	public SpectrumOperation getSpectrumOperation() {
		return OP;
	}

	
}






