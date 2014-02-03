package srp.spectrum.operator;

import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class RecombinationSpectrumOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = RecombinationSpectrumOperator.class.getSimpleName();
	public static final SpectrumOperation OP = SpectrumOperation.RECOMBINATION;
//    private Parameter parameter = null;
    private final int[] parameterWeights;
//    private double delta = 0.05;
	private int swapLength;
	
	private double autoOptimize;
//	private double scaleFactor;
	
	private int[] twoPositionIndex;
	private int[] twoSpectrumIndex;

	
	
	public RecombinationSpectrumOperator(SpectrumAlignmentModel spectrumModel) {
		super(spectrumModel, CoercionMode.COERCION_OFF);
		
		twoPositionIndex = new int[2];
		twoSpectrumIndex = new int[2];
		
//		this.delta = delta;
//		this.swapLength = swapLength;
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }
        
//		scaleFactor = (int) (spectrumLength*0.01);
//
//		if (scaleFactor <1) {
//			scaleFactor = 1;
//		}
////		convertToAutoOptimize(this.swapLength);
		
	}

	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();
		
		twoSpectrumIndex[0] = MathUtils.nextInt(spectrumCount);
//		twoSpectrumIndex[1] = -1; 
				
		do {
			twoSpectrumIndex[1] = MathUtils.nextInt(spectrumCount);
        }while (twoSpectrumIndex[0] == twoSpectrumIndex[1]);

		twoPositionIndex[0] = MathUtils.nextInt(spectrumLength);//-swapLength);
		twoPositionIndex[1]= spectrumLength;//twoPositionIndex[0] + swapLength;

		Spectrum spectrum1 = spectrumModel.getSpectrum(twoSpectrumIndex[0]);
		Spectrum spectrum2 = spectrumModel.getSpectrum(twoSpectrumIndex[1]);
		SpectraParameter spectra1;
		SpectraParameter spectra2;

//		System.err.println(swapLength);
		
		for (int i = twoPositionIndex[0]; i < twoPositionIndex[1]; i++) {
			
//			int spectrumIndex = i;
//			System.err.println(spectrumIndex +"\t"+ siteIndex);
			spectra1 = spectrum1.getSpectra(i);
			spectra2 = spectrum2.getSpectra(i);
//			System.err.println(Arrays.toString(spectra1.getParameterValues()) +"\t"+ 
//					Arrays.toString(spectra2.getParameterValues()));
			for (int j = 0; j < DIMENSION; j++) {
				double tempFreq = spectra1.getFrequency(j);
				spectra1.setParameterValueQuietly(j, spectra2.getFrequency(j));
				spectra2.setParameterValueQuietly(j, tempFreq);
				
			}
//			System.err.println("After");
//			System.err.println(Arrays.toString(spectra1.getParameterValues()) +"\t"+ 
//					Arrays.toString(spectra2.getParameterValues()));
		}
		spectrumModel.setSpectrumOperationRecord(OP, twoSpectrumIndex, twoPositionIndex);
		
		spectrumModel.endSpectrumOperation();

		return 0.0;
	}


	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

	@Override
	public double getCoercableParameter() {
	    return autoOptimize;
	}

	@Override
	public void setCoercableParameter(double autoOpt) {
		convertFromAutoOptimizeToValue(autoOpt);
	}

	private void convertFromAutoOptimizeToValue(double autoOpt) {
	    	autoOptimize = autoOpt;
			swapLength =  1 + (int) Math.exp(autoOptimize);
	    }

	private double convertToAutoOptimize(int length) {
		swapLength = length;
		checkParameterIsValid();
		autoOptimize = Math.log(swapLength - 1);
	    return autoOptimize;
	}

	private void checkParameterIsValid() {
		if (swapLength > spectrumLength){
			swapLength = spectrumLength;
		}
	}
	
	@Override
	public double getRawParameter() {
		return swapLength;
	}
//	
//
//    @Override
//	public double getCoercableParameter() {
////    	double t = Math.log(delta/(1-delta));
////    	return t;
////        return Math.log(1.0 / delta - 1.0);
//        return Math.log(delta);
//    }
//
//    @Override
//	public void setCoercableParameter(double value) {
////    	mm++;
////    	if(mm%1000 == 0){
////    		System.out.println(value +"\t"+ delta +"\t"+ getAcceptanceProbability());
////    	}
//        delta = Math.exp(value);
////        double t = Math.exp(value);
////        delta = t/(t+1);
//
//    }

//    @Override
//	public double getRawParameter() {
//        return delta;
//    }

    @Override
	public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    @Override
	public final String getPerformanceSuggestion() {
    	String s = "Tuning "+swapLength; 
    	return s;

    }

    @Override
	public String toString() {
        return getOperatorName() + "(windowsize=" + swapLength + ")";
    }


	@Override
	public SpectrumOperation getSpectrumOperation() {
		return OP;
	}
	
}



