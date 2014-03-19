package srp.spectrum.operator;

import java.util.Arrays;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.util.FastMath;

import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class SpiltFreqMultiSpectrumOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = SpiltFreqMultiSpectrumOperator.class.getSimpleName();
	public static final SpectrumOperation OP = SpectrumOperation.DELTA_MULTI;
	
    private final int[] parameterWeights;
    private double delta;
    private int swapBasesCount;
    private double autoOptimize;
    
//    private double[] debugList = new double[8];
//	private int scaleFactor=1;
//	public static double totalCount = 0;
	
	public SpiltFreqMultiSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
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

		spectrumModel.startSpectrumOperation();

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
		
		double scaleFactor = 0.9; 
		double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

		for (int i = 0; i < swapBasesCount; i++) {
//			scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
//			siteIndexs[i] = MathUtils.nextInt(spectrumLength);
//			System.err.println(spectrumIndex +"\tMultiOp\t"+ i +"\t"+ siteIndexs[i]);
			
			
			spectra[i] = spectrum.getSpectra(siteIndexs[i]);
			
				
	        
			int maxIndex = -1;
			maxIndex = MathUtils.nextInt(DIMENSION);
//			double maxFreq = 0;
			double[] oldFreq = new double[DIMENSION];
			for (int j = 0; j < DIMENSION; j++) {
				oldFreq[j] = spectra[i].getFrequency(j);
//				if(oldFreq[j] > maxFreq){
//					maxFreq = oldFreq[j];
//					maxIndex = j;
//				}
			}
			//scaleFactor large > large variance on scale
//			System.out.println(scale);
			double delta = 0.1;
//			boolean upDown = MathUtils.nextBoolean();
//			if(upDown){
//				scale = 1/scale;
//				delta = -delta;
//			}
			double delta3 = delta/3;
//			System.out.println(Arrays.toString(oldFreq) +"\t"+ scale);
			double sum = 0;
			double[] newFreq = new double[DIMENSION];
			
			BetaDistribution b = new BetaDistribution(4*oldFreq[maxIndex], 4*(1-oldFreq[maxIndex]));
//			System.out.println( b.sample() +"\t"+ oldFreq[maxIndex] );
			b.sample();
			
			
//			Math.log(1.0 / scale - 1.0);
//	    
//
//	    public void setCoercableParameter(double value) {
//	        scaleFactor = 1.0 / (Math.exp(value) + 1.0);

			for (int j = 0; j < DIMENSION; j++) {
//				oldFreq[j] *= scale;
				if(j==maxIndex){
					newFreq[j] = oldFreq[j]* scale;
////					oldFreq[j] += delta;
				}
				else{
					newFreq[j] = oldFreq[j]/ scale;
////					oldFreq[j] -= delta3;
				}
				
				boolean upDown = MathUtils.nextBoolean();
//				 d = MathUtils.nextDouble() * delta * scalar1; //TODO test this
				if(upDown){
					newFreq[j] = oldFreq[j]* scale;
					newFreq[j] = (oldFreq[j]-0.25)* scale + 0.25; // move towards 0.25
////					oldFreq[j] += delta;
				}
				else{
					newFreq[j] = oldFreq[j]/ scale;
					newFreq[j] = (oldFreq[j]-0.25)/ scale + 0.25;
////					oldFreq[j] -= delta3;
				}
				
				
//				if (oldFreq[j]< 0.001) {
////					oldFreq[j] = 0.001;
//					System.out.println("S: "+Arrays.toString(oldFreq));
//				}
//				if (oldFreq[j]> 0.997) {
////					oldFreq[j] = 0.997;
//					System.out.println("L: "+Arrays.toString(oldFreq));
//				}
				sum += newFreq[j];
			}
//			System.out.println(Arrays.toString(oldFreq));
//			System.out.println(Arrays.toString(newFreq));
//			System.out.println(sum);
			
			for (int j = 0; j < DIMENSION; j++) {
				newFreq[j] = newFreq[j]/sum;
			}
//			ScaleOperator
//			double[] temp = new double[4]; 
//			System.arraycopy(newFreq, 0, temp, 0, 4);
			boolean isPrint = false;
			
			for (int j = 0; j < DIMENSION; j++) {
				
				if (newFreq[j]< 0.001) {
//					System.out.println("S: "+Arrays.toString(newFreq));
					newFreq[j] = 0.001;
					isPrint = true;
					
				}
				if (newFreq[j]> 0.99) {
//					System.out.println("L: "+Arrays.toString(newFreq));
					newFreq[j] = 0.99;
					isPrint = true;
				}
				
			}
			sum = 0;
			for (int j = 0; j < DIMENSION; j++) {
				sum += newFreq[j];
				
			}
//			System.out.println(sum);
			for (int j = 0; j < DIMENSION; j++) {
//				newFreq[j] = oldFreq[j]/sum;
				newFreq[j] = newFreq[j]/sum;
				spectra[i].setFrequency(j, newFreq[j]);
			}
//			if(isPrint){
//			System.out.println(Arrays.toString(oldFreq));
//			System.out.println(Arrays.toString(temp));
//			System.out.println(Arrays.toString(newFreq));
//	        System.out.println();
//			}
//	
//	        if (scalar1[i] < BOUNDS_LOWER ||
//	                scalar1[i] > BOUNDS_UPPER ||
//	                scalar2[i] < BOUNDS_LOWER ||
//	                scalar2[i] > BOUNDS_UPPER ) {
////	        	System.err.println("throw");
//	            throw new OperatorFailedException("proposed values out of range!");
//	        }
//	        
		}
//		for (int i = 0; i < swapBasesCount; i++) {
//			
//			spectra[i].setFrequency(dim1[i], scalar1[i]);
//			spectra[i].setFrequency(dim2[i], scalar2[i]);
//
//		}
        // symmetrical move so return a zero hasting ratio
		spectrumModel.setSpectrumOperationRecord(OP, spectrumIndex, siteIndexs);
		
		spectrumModel.endSpectrumOperation();
//		(goingUp - goingDown - 2) * Math.log(scale);
		
		return 2 * Math.log(scale);
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





