package srp.haplotypes.likelihood;

import java.util.ArrayList;

import dr.inference.distribution.LogLinearModel;

public class LikelihoodScaler {

	/*
	 * A/C = e^log(A/C) = e^(logA - logC) = e^B
	 * 
	 * usage: 
	 * scale tiny non-log likelihood in summation to avoid == 0
	 * 
	 * log ( sum ( tiny non-log likelihood) )  
	 * scale it with LOG_SCALE, prob/LOG_SCALE
	 * log (sum (scaled likelihood)) + LOG_SCALE
	 * 
	 * example:
	 * 		for (int i = 0; i < array.length; i++) { 
	 * 			prob += avoidUnderflow(logProb, LOG_C);}
	 * 		}
	 * 		logLikelihood = (Math.log(prob) + LOG_C);
	 * 
	 * 
	 */

	private double logScaler;
	private double sumScaledLikelihood;

	public LikelihoodScaler(double logScaler){
		this.logScaler = logScaler;
		reset();
	}
	
	public void reset() {
		sumScaledLikelihood = 0;		
	}

	public double getLogLikelihood(double sumScaledLikelihood){
		
		double logLikelihood = Math.log(sumScaledLikelihood) + logScaler;
		return logLikelihood;
	}

	public double getLogLikelihood(){
		
		double logLikelihood = Math.log(sumScaledLikelihood) + logScaler;
		return logLikelihood;
	}

	public void addScaleLogProb(double logProb){
		sumScaledLikelihood += scale(logProb, logScaler);
	}

	public void minusScaleLogProb(double logProb){
		sumScaledLikelihood -= scale(logProb, logScaler);
	}

	public void addScaleLogProbMulti(double logProb, int count){
		double sum = scale(logProb, logScaler) * count;
		sumScaledLikelihood += sum;
	}


	public double scale(double logProb){
		return scale(logProb, logScaler);
	}

	private static double scale(double logProb, double logScaler) {

		double expB = Math.exp(logProb - logScaler);

		return expB;
	}


	
	public void updateScaledLogProb(double oldScaledLogProb, double newScaledLogProb) {
		sumScaledLikelihood = sumScaledLikelihood - oldScaledLogProb + newScaledLogProb;
		
	}

	public void setsumScaledLikelihood(double scaled) {
		sumScaledLikelihood = scaled;
		
	}

	public void add(double scaledLogProb){
		sumScaledLikelihood += scaledLogProb;
	}
	public void minus(double scaledLogProb) {
		sumScaledLikelihood -= scaledLogProb;
		
	}

	public double getSumScaledLikelihood() {
		return sumScaledLikelihood;
	}

	
	@Deprecated
	public double sumLogLikelihood(double[] scaledLogLikelihood) {
		double logLikelihood = 0;
		for (int i = 0; i < scaledLogLikelihood.length; i++) {
			logLikelihood += Math.log(scaledLogLikelihood[i]) ;
		}
		logLikelihood += (logScaler*scaledLogLikelihood.length);
		return logLikelihood;
	}
}
