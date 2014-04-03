package srp.likelihood.stateLikelihood;

import dr.math.distributions.BetaDistribution;

public class BetaMeanStateLikelihood extends StateLikelihood {

	private BetaDistribution betaDMean;

	public BetaMeanStateLikelihood() {
		this(NOT_ERROR_RATE, ERROR_RATE);
	}

	public BetaMeanStateLikelihood(double alpha, double beta) {
		super();
		betaDMean = new BetaDistribution(alpha, beta);
	}

	@Override
	public double caluclateStateLogLikelihood(double frequency) {
		double logLikelihood = betaDMean.logPdf(frequency);
		return logLikelihood;
	}

}
