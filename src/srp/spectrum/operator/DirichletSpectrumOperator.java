package srp.spectrum.operator;

import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;

import org.apache.commons.math3.random.RandomDataGenerator;

import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;

import com.google.common.primitives.Ints;

import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.GammaFunction;
import dr.math.MathUtils;
import dr.math.distributions.DirichletDistribution;
import dr.math.distributions.GammaDistribution;

public class DirichletSpectrumOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = DirichletSpectrumOperator.class.getSimpleName();
//	public static final SpectrumOperation OP = SpectrumOperation.DIRICHLET;
	public static final SpectrumOperation OP = SpectrumOperation.DELTA_MULTI;
	
	private static double MIN_FREQ = 0.01;
	
//    private Parameter parameter = null;

	private static final int MIN_BASE = 1;
    private final int[] parameterWeights;
//    private double delta = 0.05;
    private int baseCount = 1;
	
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
		this.baseCount = baseCount;
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }
        
        alpha = 100;
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
	private int scaleFactor=1;
	
	
	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();

//		spectrumModel.swapHaplotypeSingleBase(OP, posChar);
		int spectrumIndex = MathUtils.nextInt(spectrumCount);
//		int siteIndex = MathUtils.nextInt(spectrumLength);

        SpectraParameter[] spectra = new SpectraParameter[baseCount];
        double[] delta = new double[baseCount];
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
		Set<Integer> generated = new HashSet<Integer>();
		while (generated.size() < baseCount)
		{
		    Integer next = MathUtils.nextInt(spectrumLength);
		    generated.add(next);
		}

		siteIndexs = Ints.toArray(generated);
		double ratio = 0;
		for (int i = 0; i < baseCount; i++) {
			
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
    	baseCount =  MIN_BASE + (int) Math.exp(autoOptimize*scaleFactor);

    	checkParameterIsValid();
		
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





