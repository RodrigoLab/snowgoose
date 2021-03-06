package srp.operator.spectrum;

import org.apache.commons.math3.util.FastMath;

import srp.evolution.OperationType;
import srp.evolution.spectrum.SpectraParameter;
import srp.evolution.spectrum.Spectrum;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class DeltaMixMultiSpectrumOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = DeltaMixMultiSpectrumOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.MULTI;
	
    private final int[] parameterWeights;
    private double delta;
    private int swapBasesCount;
    private double autoOptimize;
    
//    private double[] debugList = new double[8];
//	private int scaleFactor=1;
//	public static double totalCount = 0;
	
	public DeltaMixMultiSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			double delta, int baseCount, CoercionMode mode) {
		super(spectrumModel, mode);
		
		
		this.delta = delta;
		this.swapBasesCount = baseCount;
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }
        convertToAutoOptimize(this.swapBasesCount);

	}

	
	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startAlignmentModelOperation();

//		spectrumModel.swapHaplotypeSingleBase(OP, posChar);
		int spectrumIndex = MathUtils.nextInt(spectrumCount);
//		int siteIndex = MathUtils.nextInt(spectrumLength);

        SpectraParameter[] spectra = new SpectraParameter[swapBasesCount];
//        double[] d = new double[swapBasesCount];
//        double[] scalar1 = new double[swapBasesCount];
//        double[] scalar2 = new double[swapBasesCount];
//        int[] dim1 = new int[swapBasesCount];
//		int[] dim2 = new int[swapBasesCount];
        
//		int[] siteIndexs = new int[swapBasesCount]; 
				
		Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
		
		int[] siteIndexs = generateUniqueSites(swapBasesCount);
		
		double logq = 0;
		for (int i = 0; i < swapBasesCount; i++) {
			
			spectra[i] = spectrum.getSpectra(siteIndexs[i]);
			
				
	        // get two dimensions

	        int dim1 = MathUtils.nextInt(DIMENSION);
	        int dim2 = getAnotherDimension(dim1);
    
	        double scalar1 = spectra[i].getFrequency(dim1);
	        double scalar2 = spectra[i].getFrequency(dim2);
	
//	        d[i] = MathUtils.nextDouble() * delta;
//	        d[i] = MathUtils.nextDouble() * delta * scalar1[i]; 
//	        d[i] = 0.95 * scalar1[i];
////	        d[i] = delta;
//	        scalar1[i] -= d[i];
//	        scalar2[i] += d[i];
	
//	        System.out.println(scalar1 +"\t"+ scalar2);
	
	        final double d = MathUtils.nextDouble() * delta * scalar1;
	        scalar1 -= d;
            scalar2 += d;

            spectra[i].setFrequency(dim1, scalar1);
            spectra[i].setFrequency(dim2, scalar2);

	        logq += Math.log(scalar2 / (scalar1 + d));
	        
	        if (scalar1 < BOUNDS_LOWER ||
	                scalar1 > BOUNDS_UPPER ||
	                scalar2 < BOUNDS_LOWER ||
	                scalar2 > BOUNDS_UPPER ) {
//	        	System.err.println("throw");
	        	System.out.println(scalar1 +"\t"+ scalar2);
	            throw new OperatorFailedException("proposed values out of range!");
	        }
//	        
		}
//		for (int i = 0; i < swapBasesCount; i++) {
//			
//			spectra[i].setFrequency(dim1[i], scalar1[i]);
//			spectra[i].setFrequency(dim2[i], scalar2[i]);
//
//		}
        // symmetrical move so return a zero hasting ratio
		spectrumModel.setOperationRecord(OP, spectrumIndex, siteIndexs);
		
		spectrumModel.endAlignmentModelOperation();

		return logq;
	}


	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}


//    @Override
//	public double getCoercableParameter() {
//        return Math.log(delta);
//    }
//
//    @Override
//	public void setCoercableParameter(double value) {
//        delta = Math.exp(value);
//
//    }

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
	    	swapBasesCount =  MIN_BASE + (int) FastMath.exp(autoOptimize);
//			System.out.println(autoOptimize +"\t"+ Math.exp(autoOptimize*scaleFactor));
			
//			System.out.print("A=" + swapLength + "\t" + autoOptimize + "\t" +
//					"accept: " + getAcceptCount()/(double)getCount() + "\t"  );
			
			checkParameterIsValid();
			
//			System.out.print("newL:"+swapLength+" ");
	//		System.out.print("A\t" + swapFragmentLength + "\t" + autoOptimize + "\t"  );
	    }

	private double convertToAutoOptimize(int length) {
		swapBasesCount = length;
		checkParameterIsValid();
		autoOptimize = Math.log(swapBasesCount - MIN_BASE);
	    return autoOptimize;
	}

	private void checkParameterIsValid() {
		if (swapBasesCount > spectrumLength){
			swapBasesCount = spectrumLength;
		}
	}
	
    @Override
	public double getRawParameter() {
//        return delta;
        return swapBasesCount;
    }

    @Override
	public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    @Override
	public final String getPerformanceSuggestion() {
    	String s = "Tuning delat:"+delta; 
    	return s;

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





