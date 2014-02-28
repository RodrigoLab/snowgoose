package srp.spectrum.operator;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.FastMath;

import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;

import com.google.common.primitives.Ints;
import com.sun.org.glassfish.external.statistics.Stats;

import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.GammaFunction;
import dr.math.MathUtils;

public class DirichletSpectrumOperator extends AbstractDirichletSpectrumOperator{

	public static final String OPERATOR_NAME = DirichletSpectrumOperator.class.getSimpleName();
//	public static final SpectrumOperation OP = SpectrumOperation.DIRICHLET;
	public static final SpectrumOperation OP = SpectrumOperation.DELTA_MULTI;
	
	private static final int MIN_BASE = 1;
	private static final double MIN_FREQ = SpectraParameter.MIN_FREQ;
	private static final double MAX_FREQ = SpectraParameter.MAX_FREQ;
	
    private final int[] parameterWeights;


    private int swapBasesCount;
    private double alpha;
	
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

		int spectrumIndex = MathUtils.nextInt(spectrumCount);
		Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);

		int[] siteIndexs = generateSiteIndexs(swapBasesCount, spectrumLength);
		double ratio = 0;
		for (int i = 0; i < swapBasesCount; i++) {
			
			SpectraParameter spectra = spectrum.getSpectra(siteIndexs[i]);
			double sum = 0;

			for (int j = 0; j < newFreq.length; j++) {
				oldFreq[j] = spectra.getFrequency(j);
				if(oldFreq[j]==0){
					oldFreq[j] = MIN_FREQ;
				}
				oldParameter[j] = oldFreq[j]*alpha;

				newFreq[j] = MathUtils.nextGamma(oldParameter[j], 1);
				sum += newFreq[j]; 
			}
			for (int j = 0; j < newFreq.length; j++) {
				newFreq[j] /= sum;
				if(newFreq[j]<MIN_FREQ){
					newFreq[j] = MIN_FREQ;
				}
				else if(newFreq[j]>MAX_FREQ){
					newFreq[j] = MAX_FREQ;
				}
			}
			
			for (int j = 0; j < newFreq.length; j++) {
				newParameter[j] = newFreq[j]*alpha;
				spectra.setFrequency(j, newFreq[j]);
			}
			
//			 get proposal ratio
//			DirichletUtils.calculatelogq(oldFreq, newFreq, oldParameter, newParameter);			
			double x = dirichletLnPdf(newParameter, oldFreq);
			double y = dirichletLnPdf(oldParameter, newFreq);
			ratio += (x - y);
			

		}

		spectrumModel.setSpectrumOperationRecord(OP, spectrumIndex, siteIndexs);
		
		spectrumModel.endSpectrumOperation();
//		System.out.print("diriMulti: "+ratio +"\t");
		return ratio;
	}
	/* get prior ratio 
	x = y = 0.0;
	for (i=0; i<nRates; i++)
		x += (alphaDir[i]-1.0)*log(newRate[i]);
	for (i=0; i<nRates; i++)
		y += (alphaDir[i]-1.0)*log(oldRate[i]);
	(*lnPriorRatio) = x - y;
	*/
	





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
	
	public double getAlpha(){
		return alpha;
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
