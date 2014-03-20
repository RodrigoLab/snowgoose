package srp.spectrum.operator;

import srp.spectrum.SpectraParameter;
import srp.spectrum.SpectrumAlignmentModel;
import dr.inference.operators.CoercionMode;
import dr.math.MathUtils;


public abstract class AbstractSwapSpectrumOperator extends AbstractSpectrumOperator{

	protected boolean random = true;
	public AbstractSwapSpectrumOperator(SpectrumAlignmentModel spectrumModel,
			CoercionMode mode, boolean random) {
		super(spectrumModel, mode);
		this.random = random;
	}
	/**
	 * random==true: any two dimension
	 * random==falso: and 1 only
	 * @param spectra
	 */
	
	public void swapFrequency(SpectraParameter spectra) {
		swapFrequency(spectra, random);
		
//		SpectraParameter.checkSpectra(spectra);
//        int dim1 = MathUtils.nextInt(DIMENSION);
//        int dim2;// = dim1;
//        do {
//            dim2 = MathUtils.nextInt(DIMENSION);
//        }while (dim1 == dim2);
//        
//        double scalar1 = spectra.getParameterValue(dim1);
//        double scalar2 = spectra.getParameterValue(dim2);
//        spectra.setParameterValue(dim1, scalar2);
//        spectra.setParameterValue(dim2, scalar1);
	}

	public void swapFrequency(SpectraParameter spectra, boolean random) {
		if(random){
			int dim1 = MathUtils.nextInt(DIMENSION);
			int dim2 = getAnotherDimension(dim1);
	        
	        double scalar1 = spectra.getFrequency(dim1);
	        double scalar2 = spectra.getFrequency(dim2);
	        spectra.setFrequency(dim1, scalar2);
	        spectra.setFrequency(dim2, scalar1);	
		}

		else {//get Max
//			double tempValue = -1;
//			double maxDim = -1;
			int dim1 = 0;
			double scalar1 = spectra.getParameterValue(dim1);;
			for (int d = 1; d < DIMENSION; d++) {
				double tempValue = spectra.getParameterValue(d);
				if (tempValue > scalar1) {
					dim1 = d;
					scalar1 = tempValue;
				}
			}
			int dim2 = getAnotherDimension(dim1);
			double scalar2 = spectra.getParameterValue(dim2);        
			
	        spectra.setFrequency(dim1, scalar2);
	        spectra.setFrequency(dim2, scalar1);	

			
		}
		
	}

}
