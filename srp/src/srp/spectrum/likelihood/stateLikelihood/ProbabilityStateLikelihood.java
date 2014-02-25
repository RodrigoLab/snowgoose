package srp.spectrum.likelihood.stateLikelihood;

public class ProbabilityStateLikelihood extends StateLikelihood {
	
	public ProbabilityStateLikelihood(){
		super();
	}
	
	@Override
	public double caluclateStateLogLikelihood(double frequency) {
		double logLikelihood = Math.log(frequency * NOT_ERROR_RATE
				+ (1 - frequency) * ERROR_RATE);
		return logLikelihood;
	}

}
