package srp.operator.spectrum;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.util.FastMath;

import srp.evolution.OperationType;
import srp.evolution.spectrum.SpectraParameter;
import srp.evolution.spectrum.Spectrum;
import srp.evolution.spectrum.SpectrumAlignmentModel;

import com.google.common.primitives.Ints;

import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.GammaFunction;
import dr.math.MathUtils;

public class DirichletAlphaSpectrumOperator extends AbstractDirichletSpectrumOperator {

	
	public static final String OPERATOR_NAME = DirichletAlphaSpectrumOperator.class.getSimpleName();
//	public static final SpectrumOperation OP = SpectrumOperation.DIRICHLET;
	public static final OperationType OP = OperationType.SINGLE;
	
//    private Parameter parameter = null;


    private final int[] parameterWeights;
//    private double delta = 0.05;
//    private final int swapBasesCount;
	
    private double autoOptimize;
    private double alpha = 100;
	
	double[] oldParameter = new double[DIMENSION];
	double[] newParameter = new double[DIMENSION];
	double[] oldFreq = new double[DIMENSION];
	double[] newFreq = new double[DIMENSION];

//	int[] siteIndexs; 
//	private double[] debugList = new double[8];
	
    
	public DirichletAlphaSpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			double alpha, CoercionMode mode) {
		super(spectrumModel, mode);
		
		
//		this.delta = delta;
//		baseCount = 
		this.alpha = alpha;
//		this.siteIndexs = new int[swapBasesCount];
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
	
	
//	private int scaleFactor=1;
	
//	TimeTrial: 22998342634	2299/calculation	10000000 ite.	DirichletAlphaSpectrumOperator	Operator only
//	2349/calculation
	
	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startAlignmentModelOperation();

		int spectrumIndex = MathUtils.nextInt(spectrumCount);
		Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);

		int siteIndex = MathUtils.nextInt(spectrumLength);


		SpectraParameter spectra = spectrum.getSpectra(siteIndex);
		nextDirichlet(spectra, alpha, oldFreq, oldParameter, newFreq, newParameter);
		
//			 get proposal ratio
		double x = dirichletLnPdf(newParameter, oldFreq);
		double y = dirichletLnPdf(oldParameter, newFreq);
		
		double ratio = (x - y);
		
		spectrumModel.setOperationRecord(OP, spectrumIndex, siteIndex);
		
		spectrumModel.endAlignmentModelOperation();
//		System.out.print("diriAlpha: "+ratio +"\t");
//		System.out.println(spectrumIndex +"\t"+ Arrays.toString(siteIndexs) +"\t"+ Arrays.toString(oldFreq) +"\t"+ Arrays.toString(newFreq));
		return ratio;
	}

	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}


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
//    	alpha =  1/FastMath.exp(autoOptimize);
    	alpha =  FastMath.exp(-autoOptimize);
		
    }

	private double convertToAutoOptimize() {
//		autoOptimize = Math.log(1/alpha);
		autoOptimize = -Math.log(alpha);
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
	public OperationType getSpectrumOperation() {
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
