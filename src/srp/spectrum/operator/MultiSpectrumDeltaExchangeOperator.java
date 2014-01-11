package srp.spectrum.operator;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.Set;

import com.google.common.collect.Collections2;
import com.google.common.primitives.Ints;

import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import srp.spectrum.SpectrumOperationRecord;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class MultiSpectrumDeltaExchangeOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = MultiSpectrumDeltaExchangeOperator.class.getSimpleName();
	public static final SpectrumOperation OP = SpectrumOperation.MULTI_DELTA;
//    private Parameter parameter = null;
    private final int[] parameterWeights;
    private double delta = 0.05;
    private int baseCount = 3;
    
	public MultiSpectrumDeltaExchangeOperator(SpectrumAlignmentModel spectrumModel, 
			double delta, CoercionMode mode) {
		super(spectrumModel, mode);
		
		
		this.delta = delta;
		baseCount = (int) (spectrumLength * 0.01);
		if(baseCount<2){
			baseCount = 3;
		}

        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }

	}

	private double[] debugList = new double[8];
	
	
	
	
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();

//		spectrumModel.swapHaplotypeSingleBase(OP, posChar);
		int spectrumIndex = MathUtils.nextInt(spectrumCount);
//		int siteIndex = MathUtils.nextInt(spectrumLength);

        SpectraParameter[] spectra = new SpectraParameter[baseCount];
        double[] d = new double[baseCount];
        double[] scalar1 = new double[baseCount];
        double[] scalar2 = new double[baseCount];
        int[] dim1 = new int[baseCount];
		int[] dim2 = new int[baseCount];
        
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
		Set<Integer> generated = new LinkedHashSet<Integer>();
		while (generated.size() < baseCount)
		{
		    Integer next = MathUtils.nextInt(spectrumLength);
		    // As we're adding to a set, this will automatically do a containment check
		    generated.add(next);
		}
//		generated.toArray(siteIndex);
		siteIndexs = Ints.toArray(generated);
//		System.out.println(Arrays.toString(siteIndex));
		for (int i = 0; i < baseCount; i++) {
			
//			siteIndexs[i] = MathUtils.nextInt(spectrumLength);
//			System.err.println(spectrumIndex +"\tMultiOp\t"+ i +"\t"+ siteIndexs[i]);
			
			
			spectra[i] = spectrum.getSpectra(siteIndexs[i]);
			
			
	        // get two dimensions
	        dim1[i] = MathUtils.nextInt(DIMENSION);
	        do {
	            dim2[i] = MathUtils.nextInt(DIMENSION);
	        }while (dim1[i] == dim2[i]);
    
	        scalar1[i] = spectra[i].getParameterValue(dim1[i]);
	        scalar2[i]= spectra[i].getParameterValue(dim2[i]);
	
	        d[i] = MathUtils.nextDouble() * delta;
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
		for (int i = 0; i < baseCount; i++) {
			
			spectra[i].setParameterValue(dim1[i], scalar1[i]);
			spectra[i].setParameterValue(dim2[i], scalar2[i]);

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
	public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    @Override
	public final String getPerformanceSuggestion() {
		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		String s = record.getSpectrumIndex() +"\t"+ Arrays.toString(record.getAllSiteIndexs()) +"\n";
		s+= spectrumModel.diagnostic() +"\n";
		s += Arrays.toString(debugList);
		spectrumModel.restoreModelState();
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

	
}





