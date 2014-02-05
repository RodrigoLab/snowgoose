package srp.spectrum.operator;

import java.util.Arrays;
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

public class DirichletAlphaSpectrumOperator extends AbstractDirichletSpectrumOperator {
			 
	public static final String OPERATOR_NAME = DirichletAlphaSpectrumOperator.class.getSimpleName();
//	public static final SpectrumOperation OP = SpectrumOperation.DIRICHLET;
	public static final SpectrumOperation OP = SpectrumOperation.DELTA_MULTI;
	
	private static double MIN_FREQ = 0.01;
	
//    private Parameter parameter = null;

	private static final int MIN_BASE = 1;
    private final int[] parameterWeights;
//    private double delta = 0.05;
    private final int swapBasesCount=1;
	
    private double autoOptimize;
    private double alpha = 100;
	
	double[] oldParameter = new double[DIMENSION];
	double[] newParameter = new double[DIMENSION];
	double[] oldFreq = new double[DIMENSION];
	double[] newFreq = new double[DIMENSION];

	int[] siteIndexs; 

    
	public DirichletAlphaSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			double alpha, CoercionMode mode) {
		super(spectrumModel, mode);
		
		
//		this.delta = delta;
//		baseCount = 
		this.alpha = alpha;
		this.siteIndexs = new int[swapBasesCount];
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }
        
//        alpha = 100;
        convertToAutoOptimize();
//        double[] alphas = new double[]{alpha, alpha, alpha, alpha};
//        DirichletDistribution dd = new DirichletDistribution(alphas);
//        GammaDistribution drGamma = new GammaDistribution(alpha, 1);
//        drGamma.nextGamma();
//        RandomDataGenerator rd = new RandomDataGenerator();
//        org.apache.commons.math3.distribution.GammaDistribution commonGamma = new org.apache.commons.math3.distribution.GammaDistribution(alpha, 1);
//        commonGamma.sample();
        

	}
	private double[] debugList = new double[8];
	
//	private int scaleFactor=1;
	
	
	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();

		int spectrumIndex = MathUtils.nextInt(spectrumCount);
		Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);

		siteIndexs[0] = MathUtils.nextInt(spectrumLength);
		double ratio = 0;

		SpectraParameter spectra = spectrum.getSpectra(siteIndexs[0]);
		double sum = 0;

		for (int j = 0; j < newFreq.length; j++) {
			oldFreq[j] = spectra.getFrequency(j);
			if(oldFreq[j]==0){
				oldFreq[j] = MIN_FREQ;
			}
			oldParameter[j] = oldFreq[j]*alpha;
			newFreq[j] = MathUtils.nextGamma(oldParameter[j], 1);
			if(newFreq[j]<MIN_FREQ){
				newFreq[j] = MIN_FREQ;
			}
			sum += newFreq[j]; 
		}
		for (int j = 0; j < newFreq.length; j++) {
			newFreq[j] /= sum;
			newParameter[j] = newFreq[j]*alpha;
			spectra.setParameterValue(j, newFreq[j]);
		}
		
		
//			 get proposal ratio
		double x = dirichletLnPdf(newParameter, oldFreq);
		double y = dirichletLnPdf(oldParameter, newFreq);
		
		ratio += (x - y);
		
		spectrumModel.setSpectrumOperationRecord(OP, spectrumIndex, siteIndexs);
		
		spectrumModel.endSpectrumOperation();
//		System.out.print("diriAlpha: "+ratio +"\t");
//		System.out.println(spectrumIndex +"\t"+ Arrays.toString(siteIndexs) +"\t"+ Arrays.toString(oldFreq) +"\t"+ Arrays.toString(newFreq));
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
//    	swapBasesCount =  MIN_BASE + (int) Math.exp(autoOptimize);
    	alpha =  Math.exp(autoOptimize);
//    	checkParameterIsValid();
		
    }

	private double convertToAutoOptimize() {
//		swapBasesCount = length;
//		checkParameterIsValid();
		autoOptimize = Math.log(alpha);
	    return autoOptimize;
	}

	
	
    @Override
	public double getRawParameter() {
//        return delta;
        return alpha;
    }

    @Override
	public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    @Override
	public final String getPerformanceSuggestion() {
    	String s = "Tuning alpha: "+alpha; 
    	return s;

    }

    @Override
	public String toString() {
        return getOperatorName() + "(Alpah=" + alpha + ")";
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
