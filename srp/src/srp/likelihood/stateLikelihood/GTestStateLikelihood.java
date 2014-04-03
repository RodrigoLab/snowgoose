package srp.likelihood.stateLikelihood;

import dr.math.distributions.GammaDistribution;

public class GTestStateLikelihood extends StateLikelihood {

	private GammaDistribution chisqD;

	public GTestStateLikelihood() {
		super();
		chisqD = new GammaDistribution(1.0 / 2.0, 2);
	}

	@Override
	public double caluclateStateLogLikelihood(double frequency) {
		double nf = 1 - frequency;
		double gt = 2 * (frequency * Math.log(frequency / NOT_ERROR_RATE) + nf
				* Math.log(nf / ERROR_RATE));
		double logLikelihood = chisqD.logPdf(gt);
		return logLikelihood;
	}

}
