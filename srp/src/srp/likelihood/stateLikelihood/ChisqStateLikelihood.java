package srp.likelihood.stateLikelihood;

import dr.math.distributions.GammaDistribution;

public class ChisqStateLikelihood extends StateLikelihood {

	private GammaDistribution chisqD;

	public ChisqStateLikelihood() {
		super();
		chisqD = new GammaDistribution(1.0 / 2.0, 2);
	}

	@Override
	public double caluclateStateLogLikelihood(double frequency) {
		double A = (frequency - NOT_ERROR_RATE);
		double B = (1 - frequency) - ERROR_RATE;
		double chi = (A * A / NOT_ERROR_RATE) + (B * B / ERROR_RATE);
		double logLikelihood = chisqD.logPdf(chi);
		return logLikelihood;
	}

}
