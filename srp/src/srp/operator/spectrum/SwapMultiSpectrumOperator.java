package srp.operator.spectrum;

import srp.evolution.OperationType;
import srp.evolution.spectrum.SpectraParameter;
import srp.evolution.spectrum.Spectrum;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class SwapMultiSpectrumOperator extends AbstractSwapSpectrumOperator {

	public static final String OPERATOR_NAME = SwapMultiSpectrumOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.MULTI;

    private final int[] parameterWeights;
//    private double delta = 0.05;
    
    private int swapBasesCount;
    private double autoOptimize;
    
	public SwapMultiSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			int baseCount, CoercionMode mode) {
		this(spectrumModel, baseCount, true, mode);
	}
	public SwapMultiSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			int baseCount, boolean random, CoercionMode mode) {
		super(spectrumModel, mode, random);
		
		
//		this.delta = delta;
		this.swapBasesCount = baseCount;
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }
        convertToAutoOptimize(this.swapBasesCount);

	}

//	private double[] debugList = new double[8];
//	private int scaleFactor=1;
//	public static int totalCount = 0;

	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startAlignmentModelOperation();

//		spectrumModel.swapHaplotypeSingleBase(OP, posChar);
		int spectrumIndex = MathUtils.nextInt(spectrumCount);
//		int siteIndex = MathUtils.nextInt(spectrumLength);

//        SpectraParameter[] spectra = new SpectraParameter[baseCount];
//        double[] d = new double[baseCount];
//        double[] scalar1 = new double[baseCount];
//        double[] scalar2 = new double[baseCount];
//        int[] dim1 = new int[baseCount];
//		int[] dim2 = new int[baseCount];
        
//		int[] siteIndexs = new int[swapBasesCount]; 
				
		Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
		
//		List<Integer> list=new ArrayList<Integer>();
//	    while(count<50){
//	        int num=random.nextInt(50);
//	            if(!list.contains(num)){
//	                list.add(num);
//	                ++count;  
//	            }                    
//	    }
//	    
//	    
		int[] siteIndexs = 
				generateUniqueSites(swapBasesCount);
//		int[] siteIndexs = randomSampleSites(swapBasesCount);
		//		System.out.println(Arrays.toString(siteIndex));
		for (int i = 0; i < swapBasesCount; i++) {
			
			SpectraParameter spectra = spectrum.getSpectra(siteIndexs[i]);
			swapFrequency(spectra);
	        
		}
        // symmetrical move so return a zero hasting ratio
		spectrumModel.setOperationRecord(OP, spectrumIndex, siteIndexs);
		
		spectrumModel.endAlignmentModelOperation();

		return 0.0;
	}

	
	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME+"(random="+random+")";
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
//	    	baseCount =  MIN_BASE + (int) Math.exp(autoOptimize*scaleFactor);
	    	swapBasesCount =  MIN_BASE + (int) Math.pow(2, autoOptimize);
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
//		autoOptimize = Math.log(baseCount - MIN_BASE)/scaleFactor;
		autoOptimize = Math.sqrt(swapBasesCount - MIN_BASE);
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
    	String s = "Tuning "+swapBasesCount; 
    	return s;

    }

    @Override
	public String toString() {
        return getOperatorName() + "(windowsize=" + swapBasesCount + ")";
    }


	@Override
	public OperationType getSpectrumOperation() {
		return OP;
	}

	
}





