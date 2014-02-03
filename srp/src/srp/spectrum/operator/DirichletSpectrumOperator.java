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
import dr.math.GammaFunction;
import dr.math.MathUtils;

public class DirichletSpectrumOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = DirichletSpectrumOperator.class.getSimpleName();
//	public static final SpectrumOperation OP = SpectrumOperation.DIRICHLET;
	public static final SpectrumOperation OP = SpectrumOperation.DELTA_MULTI;
	
	private static double MIN_FREQ = 0.01;
	
//    private Parameter parameter = null;

	private static final int MIN_BASE = 1;
    private final int[] parameterWeights;
//    private double delta = 0.05;
    private int swapBasesCount;
	
    private double alpha = 100;
	
	double[] oldParameter = new double[DIMENSION];
	double[] newParameter = new double[DIMENSION];
	double[] oldFreq = new double[DIMENSION];
	double[] newFreq = new double[DIMENSION];

    
	public DirichletSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			int baseCount, CoercionMode mode) {
		super(spectrumModel, mode);
		
		
//		this.delta = delta;
//		baseCount = 
		this.swapBasesCount = baseCount;
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }
        
        alpha = 100;
        convertToAutoOptimize(this.swapBasesCount);
//        double[] alphas = new double[]{alpha, alpha, alpha, alpha};
//        DirichletDistribution dd = new DirichletDistribution(alphas);
//        GammaDistribution drGamma = new GammaDistribution(alpha, 1);
//        drGamma.nextGamma();
//        RandomDataGenerator rd = new RandomDataGenerator();
//        org.apache.commons.math3.distribution.GammaDistribution commonGamma = new org.apache.commons.math3.distribution.GammaDistribution(alpha, 1);
//        commonGamma.sample();
        

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
        double[] delta = new double[swapBasesCount];
        double[] scalar1 = new double[swapBasesCount];
        double[] scalar2 = new double[swapBasesCount];
        int[] dim1 = new int[swapBasesCount];
		int[] dim2 = new int[swapBasesCount];
        
		int[] siteIndexs = new int[swapBasesCount]; 
				
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
		while (generated.size() < swapBasesCount)
		{
		    Integer next = MathUtils.nextInt(spectrumLength);
		    generated.add(next);
		}

		siteIndexs = Ints.toArray(generated);
		double ratio = 0;
		for (int i = 0; i < swapBasesCount; i++) {
			
//			siteIndexs[i] = MathUtils.nextInt(spectrumLength);
//			System.err.println(spectrumIndex +"\tMultiOp\t"+ i +"\t"+ siteIndexs[i]);
			spectra[i] = spectrum.getSpectra(siteIndexs[i]);
			double sum = 0;

//			System.out.println(Arrays.toString(spectra[i].getFrequencies()));
			boolean isNotValid = true;
			for (int j = 0; j < oldParameter.length; j++) {
				oldFreq[j] = spectra[i].getFrequency(j);
				oldParameter[j] = oldFreq[j]*alpha;
			}
			do{
				isNotValid = false;
				sum = 0;
				
				for (int j = 0; j < newFreq.length; j++) {
					//TODO remove later
					if(oldFreq[j]< MIN_FREQ){
						oldFreq[j] = MIN_FREQ;
						oldParameter[j] = oldFreq[j]*alpha;
//						oldParameter[j]=0.1;
					}
					newFreq[j] = MathUtils.nextGamma(oldParameter[j], 1);
					sum += newFreq[j]; 
				}
				for (int j = 0; j < newFreq.length; j++) {
					newFreq[j] /= sum;
					if(newFreq[j]< MIN_FREQ ){
						isNotValid = true;
						break;
					}
					
				}
			}while(isNotValid);
			for (int j = 0; j < newFreq.length; j++) {
				spectra[i].setParameterValue(j, newFreq[j]);
			}
			
			
			

//			 get proposal ratio 
			sum = 0.0;
			for (int d=0; d<DIMENSION; d++){
				newParameter[d] = newFreq[d]*alpha;
				sum += newParameter[d];
			}
			double x = GammaFunction.lnGamma(sum);
			for (int d=0; d<DIMENSION; d++)
				x -= GammaFunction.lnGamma(newParameter[d]);
			for (int d=0; d<DIMENSION; d++)
				x += (newParameter[d]-1.0)*Math.log(oldFreq[d]);

			sum = 0.0;
			for (int d=0; d<DIMENSION; d++)
				sum += oldParameter[d];
			double y = GammaFunction.lnGamma(sum);
			for (int d=0; d<DIMENSION; d++)
				y -= GammaFunction.lnGamma(oldParameter[d]);
			for (int d=0; d<DIMENSION; d++)
				y += (oldParameter[d]-1.0)*Math.log(newFreq[d]);

			ratio += (x - y);
		}
//		for (int i = 0; i < baseCount; i++) {
//			
//			spectra[i].setParameterValue(dim1[i], scalar1[i]);
//			spectra[i].setParameterValue(dim2[i], scalar2[i]);
//
//		}


	
//	System.out.println(x +"\t"+ y +"\t"+ ratio);		
//    GammaDistribution drGamma = new GammaDistribution(alpha, 1);


//	for (int i=0; i<DIMENSION; i++){
////		drGamma.setShape(oldParameter[i]);
////		sum += drGamma.logPdf(newParameter[i]);
//	    DirichletDistribution dd = new DirichletDistribution(new double[]{oldParameter[i]});
//	    sum += dd.logPdf(new double[]{newParameter[i]});
//	}
//	System.out.println(sum);
	/* get prior ratio 
	x = y = 0.0;
	for (i=0; i<nRates; i++)
		x += (alphaDir[i]-1.0)*log(newRate[i]);
	for (i=0; i<nRates; i++)
		y += (alphaDir[i]-1.0)*log(oldRate[i]);
	(*lnPriorRatio) = x - y;
	*/
		
		spectrumModel.setSpectrumOperationRecord(OP, spectrumIndex, siteIndexs, delta);
		
		spectrumModel.endSpectrumOperation();

		return ratio;
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
    	swapBasesCount =  MIN_BASE + (int) Math.exp(autoOptimize);

    	checkParameterIsValid();
		
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
    	String s = "Tuning "+swapBasesCount; 
    	return s;

    }

    @Override
	public String toString() {
        return getOperatorName() + "(windowsize=" + swapBasesCount + ")";
    }


	@Override
	public SpectrumOperation getSpectrumOperation() {
		return OP;
	}

	
}





/**
 *
    public double getCoercableParameter() {
//	     return Math.log(1.0 / scaleFactor - 1.0);
        return Math.log(scaleFactor);
    }

    public void setCoercableParameter(double value) {
//	     scaleFactor = 1.0 / (Math.exp(value) + 1.0);
        scaleFactor = Math.exp(value);
    }
    ===================
    
    public double getCoercableParameter() {
        return Math.log(delta)-Math.log(1-delta);
    }

    public void setCoercableParameter(double value) {
        delta = 1/(1+Math.exp(-value));
    }

    ===================
    
    private int convertFromAutoOptimizeValue(double value) {
        return 1 + (int) Math.exp(autoOptimize);
    }

    private double convertToAutoOptimizeValue(int value) {
        return Math.log(value - 1);
    }

    public double getCoercableParameter() {
        return autoOptimize;
    }

    public void setCoercableParameter(double value) {
        autoOptimize = value;
    }


============================
    public double getCoercableParameter() {
//        return Math.log(scaleFactor);
        return Math.sqrt(scaleFactor - 1);
    }

    public void setCoercableParameter(double value) {
//        scaleFactor = Math.exp(value);
        scaleFactor = 1 + value * value;
    }
    ===================
    
        
            public double getCoercableParameter() {
        return Math.log(1.0 / scaleFactor - 1.0);
    }

    public void setCoercableParameter(double value) {
        scaleFactor = 1.0 / (Math.exp(value) + 1.0);
    }
    
    ==============================
    
    
    public double getCoercableParameter() {
        return Math.log(windowSize);
    }

    public void setCoercableParameter(double value) {
        windowSize = (int) Math.exp(value);
    }



=========================
    
    
    *
    */
