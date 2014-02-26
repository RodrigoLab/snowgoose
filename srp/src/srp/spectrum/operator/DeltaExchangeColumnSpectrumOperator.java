package srp.spectrum.operator;

import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class DeltaExchangeColumnSpectrumOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = DeltaExchangeColumnSpectrumOperator.class.getSimpleName();
	public static final SpectrumOperation OP = SpectrumOperation.DELTA_COLUMN;
//    private Parameter parameter = null;
    private final int[] parameterWeights;
    private double delta = 0.05;
    
	public DeltaExchangeColumnSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
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
	
	
	
	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();

//		int[] posChar = alignmentMapping.getNextBase();
//		spectrumModel.swapHaplotypeSingleBase(OP, posChar);
//		int spectrumIndex = MathUtils.nextInt(spectrumCount);
		int siteIndex = MathUtils.nextInt(spectrumLength);

    	
        SpectraParameter[] spectra = new SpectraParameter[spectrumCount];
        double[] d = new double[spectrumCount];
        double[] scalar1 = new double[spectrumCount];
        double[] scalar2 = new double[spectrumCount];
        int[] dim1 = new int[spectrumCount];
		int[] dim2 = new int[spectrumCount];
		
        
		for (int i = 0; i < spectrumCount; i++) {
			
			
			int spectrumIndex = i;
//			System.err.println(spectrumIndex +"\t"+ siteIndex);
			
			Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
			spectra[i] = spectrum.getSpectra(siteIndex);
	        // get two dimensions
	        dim1[i] = MathUtils.nextInt(DIMENSION);
	        do {
	            dim2[i] = MathUtils.nextInt(DIMENSION);
	        }while (dim1[i] == dim2[i]);
    
	        scalar1[i] = spectra[i].getFrequency(dim1[i]);
	        scalar2[i]= spectra[i].getFrequency(dim2[i]);
	
	        d[i] = MathUtils.nextDouble() * delta;
	        d[i] = delta;
	        scalar1[i] -= d[i];
	        scalar2[i] += d[i];
	
	        
	
	        if (scalar1[i] < BOUNDS_LOWER ||
	                scalar1[i] > BOUNDS_UPPER ||
	                scalar2[i] < BOUNDS_LOWER ||
	                scalar2[i] > BOUNDS_UPPER ) {
//	        	System.err.println("throw");
	            throw new OperatorFailedException("proposed values out of range!");
	        }
	        
		}
		for (int i = 0; i < spectrumCount; i++) {
			
			spectra[i].setFrequency(dim1[i], scalar1[i]);
			spectra[i].setFrequency(dim2[i], scalar2[i]);

		}
        // symmetrical move so return a zero hasting ratio
		spectrumModel.setSpectrumOperationRecord(OP, siteIndex, d);
		
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
	public final String getPerformanceSuggestion() {

    	String s = "Tuning "+delta; 
    	return s;
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
