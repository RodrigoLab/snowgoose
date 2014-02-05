package srp.spectrum.operator;

import srp.spectrum.SpectraParameter;
import srp.spectrum.SpectrumAlignmentModel;
import dr.inference.operators.CoercionMode;
import dr.math.MathUtils;


public abstract class AbstractSwapSpectrumOperator extends AbstractSpectrumOperator{

	private boolean random = true;
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
			// get any two dims and swap
	        int dim1 = MathUtils.nextInt(DIMENSION);
	        int dim2;// = dim1;
	        do {
	            dim2 = MathUtils.nextInt(DIMENSION);
	        }while (dim1 == dim2);
	        
	        double scalar1 = spectra.getParameterValue(dim1);
	        double scalar2 = spectra.getParameterValue(dim2);
	        spectra.setParameterValue(dim1, scalar2);
	        spectra.setParameterValue(dim2, scalar1);	
		}

		else {//0 and 1 only
			int dim1 = 0;
			double scalar1 = 0;
			for (int d = 0; d < DIMENSION; d++) {
				scalar1 = spectra.getParameterValue(d);
				if (scalar1 == 1) {
					dim1 = d;
					break;
				}
			}
			int dim2;// = dim1;
			do {
				dim2 = MathUtils.nextInt(DIMENSION);
			} while (dim1 == dim2);
			double scalar2 = spectra.getParameterValue(dim2);        

	        spectra.setParameterValue(dim1, scalar2);
	        spectra.setParameterValue(dim2, scalar1);	

			
		}
		
	}


}
