package srp.spectrum.operator;

import srp.spectrum.SpectrumAlignmentModel;
import dr.inference.operators.CoercionMode;
import dr.math.GammaFunction;

public abstract class AbstractDirichletSpectrumOperator extends AbstractSpectrumOperator{


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

}
