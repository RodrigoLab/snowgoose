package srp.spectrum.operator;

import java.util.HashSet;
import java.util.Set;

import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;

import com.google.common.primitives.Ints;

import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class SwapMultiSpectrumOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = SwapMultiSpectrumOperator.class.getSimpleName();
	public static final SpectrumOperation OP = SpectrumOperation.SWAP_MULTI;
	private static final int MIN_BASE = 3;
    private final int[] parameterWeights;
//    private double delta = 0.05;
    private int baseCount = 3;
    
	public SwapMultiSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			int baseCount, CoercionMode mode) {
		super(spectrumModel, mode);
		
		
//		this.delta = delta;
		this.baseCount = baseCount;
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }

	}

	private double[] debugList = new double[8];
	private double autoOptimize;
	private int scaleFactor=1;
	
	
	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();

//		spectrumModel.swapHaplotypeSingleBase(OP, posChar);
		int spectrumIndex = MathUtils.nextInt(spectrumCount);
//		int siteIndex = MathUtils.nextInt(spectrumLength);

//        SpectraParameter[] spectra = new SpectraParameter[baseCount];
//        double[] d = new double[baseCount];
//        double[] scalar1 = new double[baseCount];
//        double[] scalar2 = new double[baseCount];
//        int[] dim1 = new int[baseCount];
//		int[] dim2 = new int[baseCount];
        
		int[] siteIndexs = new int[baseCount]; 
				
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
		Set<Integer> generated = new HashSet<Integer>();
		while (generated.size() < baseCount)
		{
		    Integer next = MathUtils.nextInt(spectrumLength);
		    generated.add(next);
		}
//		generated.toArray(siteIndex);
		siteIndexs = Ints.toArray(generated);
		
		//		System.out.println(Arrays.toString(siteIndex));
		for (int i = 0; i < baseCount; i++) {
			
//			siteIndexs[i] = MathUtils.nextInt(spectrumLength);
//			System.err.println(spectrumIndex +"\tMultiOp\t"+ i +"\t"+ siteIndexs[i]);
			
			
			SpectraParameter spectra = spectrum.getSpectra(siteIndexs[i]);
			int dim1 = -1;
	        double scalar1= -1;
	        for (int d = 0; d < DIMENSION; d++) {
	        	scalar1 = spectra.getParameterValue(d);
	        	if(scalar1==1){
	        		dim1=d;
	        		break;
	        	}
			}
	        int dim2;// = dim1;
	        do {
	            dim2 = MathUtils.nextInt(DIMENSION);
	        }while (dim1 == dim2);
	        double scalar2 = spectra.getParameterValue(dim2);
	        
	    	
			spectra.setParameterValue(dim1, scalar2);
			spectra.setParameterValue(dim2, scalar1);

//	        // get two dimensions
//	        dim1[i] = MathUtils.nextInt(DIMENSION);
//	        do {
//	            dim2[i] = MathUtils.nextInt(DIMENSION);
//	        }while (dim1[i] == dim2[i]);
//    
//	        scalar1[i] = spectra[i].getParameterValue(dim1[i]);
//	        scalar2[i]= spectra[i].getParameterValue(dim2[i]);
//	
//	        d[i] = MathUtils.nextDouble() * delta;
//	        d[i] = delta;
//	        scalar1[i] -= d[i];
//	        scalar2[i] += d[i];
//	
//	
//	        if (scalar1[i] < BOUNDS_LOWER ||
//	                scalar1[i] > BOUNDS_UPPER ||
//	                scalar2[i] < BOUNDS_LOWER ||
//	                scalar2[i] > BOUNDS_UPPER ) {
////	        	System.err.println("throw");
//	            throw new OperatorFailedException("proposed values out of range!");
//	        }
	        
		}
//		for (int i = 0; i < baseCount; i++) {
//			
//			spectra[i].setParameterValue(dim1[i], scalar1[i]);
//			spectra[i].setParameterValue(dim2[i], scalar2[i]);
//
//		}
        // symmetrical move so return a zero hasting ratio
		spectrumModel.setSpectrumOperationRecord(OP, spectrumIndex, siteIndexs);
		
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
	    	baseCount =  MIN_BASE + (int) Math.exp(autoOptimize*scaleFactor);
//			System.out.println(autoOptimize +"\t"+ Math.exp(autoOptimize*scaleFactor));
			
//			System.out.print("A=" + swapLength + "\t" + autoOptimize + "\t" +
//					"accept: " + getAcceptCount()/(double)getCount() + "\t"  );
			
			checkParameterIsValid();
			
//			System.out.print("newL:"+swapLength+" ");
	//		System.out.print("A\t" + swapFragmentLength + "\t" + autoOptimize + "\t"  );
	    }

	private double convertToAutoOptimize(int length) {
		baseCount = length;
		checkParameterIsValid();
		autoOptimize = Math.log(baseCount - MIN_BASE)/scaleFactor;
	    return autoOptimize;
	}

	private void checkParameterIsValid() {
		if (baseCount > spectrumLength){
			baseCount = spectrumLength;
		}
	}
	
    @Override
	public double getRawParameter() {
//        return delta;
        return baseCount;
    }

    @Override
	public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    @Override
	public final String getPerformanceSuggestion() {
    	String s = "Tuning "+baseCount; 
    	return s;

    }

    @Override
	public String toString() {
        return getOperatorName() + "(windowsize=" + baseCount + ")";
    }


	@Override
	public SpectrumOperation getSpectrumOperation() {
		return OP;
	}

	
}





