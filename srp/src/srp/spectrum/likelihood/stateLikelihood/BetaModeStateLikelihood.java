package srp.spectrum.likelihood.stateLikelihood;

import dr.math.distributions.BetaDistribution;

public class BetaModeStateLikelihood extends StateLikelihood {

	private BetaDistribution betaDMode;

	public BetaModeStateLikelihood() {
		this(1+NOT_ERROR_RATE, 1+ERROR_RATE);
		//  alpha = 1.9893, beta = 1.0107
	}

	public BetaModeStateLikelihood(double alpha, double beta) {
		super();
		betaDMode = new BetaDistribution(alpha, beta);
	}

	@Override
	public double caluclateStateLogLikelihood(double frequency) {
		double logLikelihood = betaDMode.logPdf(frequency);
		return logLikelihood;
	}

}
