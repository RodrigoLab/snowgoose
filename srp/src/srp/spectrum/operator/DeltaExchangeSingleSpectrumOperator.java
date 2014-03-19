package srp.spectrum.operator;

import org.apache.commons.math3.util.FastMath;

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
    
    private double delta;
    
//  private int[] siteIndex;
//    private double[] debugList = new double[8];
//	private int mm;
//	public int failcount = 0;
	
	public DeltaExchangeSingleSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			double delta, CoercionMode mode) {
		super(spectrumModel, mode);
		
		
		this.delta = delta;
//		this.siteIndex = new int[1];
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
		int siteIndex = MathUtils.nextInt(spectrumLength);
		
		Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
		SpectraParameter spectra = spectrum.getSpectra(siteIndex);
        // get two dimensions

        int dim1 = MathUtils.nextInt(DIMENSION);
        int dim2 = getAnotherDimension(dim1);

        double scalar1 = spectra.getFrequency(dim1);
        double scalar2 = spectra.getFrequency(dim2);

		double d = delta;
//        double d = MathUtils.nextDouble() * delta;
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

        spectra.setFrequency(dim1, scalar1);
        spectra.setFrequency(dim2, scalar2);

        // symmetrical move so return a zero hasting ratio
		spectrumModel.setSpectrumOperationRecord(OP, spectrumIndex, siteIndex, d);
		
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
	public SpectrumOperation getSpectrumOperation() {
		return OP;
	}

	
}





