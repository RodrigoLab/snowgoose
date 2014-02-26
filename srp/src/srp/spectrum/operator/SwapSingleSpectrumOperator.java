package srp.spectrum.operator;

import java.util.Arrays;

import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class SwapSingleSpectrumOperator extends AbstractSwapSpectrumOperator{

	public static final String OPERATOR_NAME = SwapSingleSpectrumOperator.class.getSimpleName();
	public static final SpectrumOperation OP = SpectrumOperation.DELTA_SINGLE;
//    private Parameter parameter = null;
    private final int[] parameterWeights;
    private int[] siteIndexs;
    private double delta = Double.NaN;

    public SwapSingleSpectrumOperator(SpectrumAlignmentModel spectrumModel) {
		this(spectrumModel, true);
    }
	public SwapSingleSpectrumOperator(SpectrumAlignmentModel spectrumModel, boolean random) {
		super(spectrumModel, CoercionMode.COERCION_OFF, random);
		
		
//		this.delta = delta;
		this.siteIndexs = new int[1];
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }

	}

	private double[] debugList = new double[8];

	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();

		int spectrumIndex = MathUtils.nextInt(spectrumCount);
		siteIndexs[0] = MathUtils.nextInt(spectrumLength);
		
		Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
		SpectraParameter spectra = spectrum.getSpectra(siteIndexs[0]);

		swapFrequency(spectra);
        // symmetrical move so return a zero hasting ratio

		spectrumModel.setSpectrumOperationRecord(OP, spectrumIndex, siteIndexs);
		
		spectrumModel.endSpectrumOperation();

		return 0.0;
	}

	
	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}


    @Override
	public double getCoercableParameter() {
        return Double.NaN;
    }

    @Override
	public void setCoercableParameter(double value) {
    }

    @Override
	public double getRawParameter() {
        return delta;
    }


    @Override
	public final String getPerformanceSuggestion() {

    	String s = "No tuning"; 
    	return s;
    	
    }

    @Override
	public String toString() {
        return getOperatorName() + "(Single)";
    }

	@Override
	public SpectrumOperation getSpectrumOperation() {
		return OP;
	}

	
}






