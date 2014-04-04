package srp.operator.spectrum;

import srp.evolution.OperationType;
import srp.evolution.spectrum.SpectraParameter;
import srp.evolution.spectrum.Spectrum;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class SwapSingleSpectrumOperator extends AbstractSwapSpectrumOperator{

	public static final String OPERATOR_NAME = SwapSingleSpectrumOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.SINGLE;
//    private Parameter parameter = null;
    private final int[] parameterWeights;
//    private int[] siteIndexs;
    

    public SwapSingleSpectrumOperator(SpectrumAlignmentModel spectrumModel) {
		this(spectrumModel, true);
    }
	public SwapSingleSpectrumOperator(SpectrumAlignmentModel spectrumModel, boolean random) {
		super(spectrumModel, CoercionMode.COERCION_OFF, random);
		
		
//		this.delta = delta;
//		this.siteIndexs = new int[1];
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }

	}


	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();

		int spectrumIndex = MathUtils.nextInt(spectrumCount);
		int siteIndex = MathUtils.nextInt(spectrumLength);
		
		Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
		SpectraParameter spectra = spectrum.getSpectra(siteIndex);

		swapFrequency(spectra);
        // symmetrical move so return a zero hasting ratio

		spectrumModel.setOperationRecord(OP, spectrumIndex, siteIndex);
		
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
        return Double.NaN;
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
	public OperationType getSpectrumOperation() {
		return OP;
	}

	
}






