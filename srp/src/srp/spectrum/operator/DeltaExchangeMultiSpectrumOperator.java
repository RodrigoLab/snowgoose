package srp.spectrum.operator;

import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.util.FastMath;

import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;

import com.google.common.primitives.Ints;

import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class DeltaExchangeMultiSpectrumOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = DeltaExchangeMultiSpectrumOperator.class.getSimpleName();
	public static final SpectrumOperation OP = SpectrumOperation.DELTA_MULTI;
	private static final int MIN_BASE = 3;
    private final int[] parameterWeights;
    private double delta = 0.05;
    private int swapBasesCount = 3;
    
	public DeltaExchangeMultiSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
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

	private double[] debugList = new double[8];
	private double autoOptimize;
//	private int scaleFactor=1;
	
	
	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();

//		spectrumModel.swapHaplotypeSingleBase(OP, posChar);
		int spectrumIndex = MathUtils.nextInt(spectrumCount);
//		int siteIndex = MathUtils.nextInt(spectrumLength);

        SpectraParameter[] spectra = new SpectraParameter[swapBasesCount];
        double[] d = new double[swapBasesCount];
        double[] scalar1 = new double[swapBasesCount];
        double[] scalar2 = new double[swapBasesCount];
        int[] dim1 = new int[swapBasesCount];
		int[] dim2 = new int[swapBasesCount];
        
		int[] siteIndexs = new int[swapBasesCount]; 
				
		Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
		

		Set<Integer> generated = new HashSet<Integer>();
		while (generated.size() < swapBasesCount)
		{
		    Integer next = MathUtils.nextInt(spectrumLength);
		    generated.add(next);
		}
//		generated.toArray(siteIndex);
		siteIndexs = Ints.toArray(generated);
//		System.out.println(Arrays.toString(siteIndex));
		for (int i = 0; i < swapBasesCount; i++) {
			
//			siteIndexs[i] = MathUtils.nextInt(spectrumLength);
//			System.err.println(spectrumIndex +"\tMultiOp\t"+ i +"\t"+ siteIndexs[i]);
			
			
			spectra[i] = spectrum.getSpectra(siteIndexs[i]);
			
			
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
		for (int i = 0; i < swapBasesCount; i++) {
			
			spectra[i].setFrequency(dim1[i], scalar1[i]);
			spectra[i].setFrequency(dim2[i], scalar2[i]);

		}
        // symmetrical move so return a zero hasting ratio
		spectrumModel.setSpectrumOperationRecord(OP, spectrumIndex, siteIndexs, d);
		
		spectrumModel.endSpectrumOperation();

		return 0.0;
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





