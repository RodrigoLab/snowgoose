package srp.spectrum.operator;

import java.util.Arrays;

import srp.spectrum.SpectraParameter;
import srp.spectrum.SpectrumAlignmentModel;
import dr.inference.operators.CoercionMode;
import dr.math.GammaFunction;
import dr.math.MathUtils;

public abstract class AbstractDirichletSpectrumOperator extends AbstractSpectrumOperator{


	protected static final double MIN_FREQ = SpectraParameter.MIN_FREQ;
	protected static final double MAX_FREQ = SpectraParameter.MAX_FREQ;

	
	public AbstractDirichletSpectrumOperator(
			SpectrumAlignmentModel spectrumModel, CoercionMode mode) {
		super(spectrumModel, mode);
	}
	
	@Deprecated
	public static double calculatelogq(double[] oldFreq, double[] newFreq,
			double[] oldParameter, double[] newParameter) {
		
		double sum = 0.0;
		for (int d=0; d<DIMENSION; d++){
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

		double ratio = (x-y);
		return ratio;
	}

	public static double dirichletLnPdf(double[] parameter, double[] var) {
		double sum = 0.0;
		for (int d=0; d<DIMENSION ; d++){
			sum += parameter[d];
		}
		double x = GammaFunction.lnGamma(sum);
		for (int d=0; d<DIMENSION; d++)
			x -= GammaFunction.lnGamma(parameter[d]);
		for (int d=0; d<DIMENSION; d++)
			x += (parameter[d]-1.0)*Math.log(var[d]);

		return x;
	}

	/*
	 * small alpha -> less similar then oldFreq
	 */
	
	protected void nextDirichlet(SpectraParameter spectra, double alpha,
			double[] oldFreq, double[] oldParameter, double[] newFreq, double[] newParameter) {

		double sum = 0;
		for (int j = 0; j < newFreq.length; j++) {
			oldFreq[j] = spectra.getFrequency(j);
			if(oldFreq[j]<MIN_FREQ){
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

	}

}
