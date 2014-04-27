package srp.operator.spectrum;

import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

import srp.evolution.OperationType;
import srp.evolution.spectrum.SpectraParameter;
import srp.evolution.spectrum.Spectrum;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class testCategoriesChangeSpectrumOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = testCategoriesChangeSpectrumOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.MULTI;
	
    private final int[] parameterWeights;
    private double delta;
    private int swapBasesCount;
    private double autoOptimize;
    
//    private double[] debugList = new double[8];
//	private int scaleFactor=1;
//	public static double totalCount = 0;
	
	public testCategoriesChangeSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
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

	public static final double MAX = 0.97;
	public static final double HALF = 0.49;
	public static final double THIRD = 0.33;
	public static final double MIN = 0.01;
	public static double[][] twoCategories = new double[][]{
		{0.49, 0.49, 0.01, 0.01},
		{0.49, 0.01, 0.49, 0.01},
		{0.49, 0.01, 0.01, 0.49},
		{0.01, 0.49, 0.49, 0.01},
		{0.01, 0.49, 0.01, 0.49},
		{0.01, 0.01, 0.49, 0.49},
	};
	public static double[][] oneCategories = new double[][]{
		{0.97, 0.01, 0.01, 0.01},
		{0.01, 0.97, 0.01, 0.01},
		{0.01, 0.01, 0.97, 0.01},
		{0.01, 0.01, 0.01, 0.97},
	};
	public static final int TOTAL_TWO_CAT;
	public static final int TOTAL_ONE_CAT;
	static {
		TOTAL_TWO_CAT = twoCategories.length;
		TOTAL_ONE_CAT = oneCategories.length;
	}
	
	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startAlignmentModelOperation();

//		spectrumModel.swapHaplotypeSingleBase(OP, posChar);
		int spectrumIndex = MathUtils.nextInt(spectrumCount);
        SpectraParameter[] spectra = new SpectraParameter[swapBasesCount];

        Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
		int[] siteIndexs = generateUniqueSites(swapBasesCount);
		

		for (int i = 0; i < swapBasesCount; i++) {
			
			spectra[i] = spectrum.getSpectra(siteIndexs[i]);
			
			int maxIndex = -1;
			int[] equalIndexs = new int[2];
			int noCat = 0;
			double[] oldFreq = new double[DIMENSION];
			for (int j = 0; j < DIMENSION; j++) {
				oldFreq[j] = spectra[i].getFrequency(j);
				if(oldFreq[j] == MAX ){
					noCat = 1;
				}
				if(oldFreq[j] == HALF ){
					noCat = 2;
				}
			}
			
			if(noCat ==1){ // 0.97
				int index = MathUtils.nextInt(TOTAL_ONE_CAT);
				for (int j = 0; j < DIMENSION; j++) {
					spectra[i].setFrequency(j, oneCategories[index][j]);
				}
			}
			
			else if(noCat == 2){//0.49
				int index = MathUtils.nextInt(TOTAL_TWO_CAT);
				for (int j = 0; j < DIMENSION; j++) {
					spectra[i].setFrequency(j, twoCategories[index][j]);	
				}
				
			}
			else{
				throw new IllegalArgumentException(
						"something wrong with the spectrum\t" + noCat + "\t"
								+ Arrays.toString(spectra[i].getFrequencies()));
			}
			
//			SpectraParameter.checkSpectra(spectra[i]);//REMOVE
		}

		spectrumModel.setOperationRecord(OP, spectrumIndex, siteIndexs);
		
		spectrumModel.endAlignmentModelOperation();
		
		return 0;
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
	public OperationType getSpectrumOperation() {
		return OP;
	}

	
}





