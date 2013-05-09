package srp.likelihood;

import java.util.ArrayList;

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

	public double getLogLikelihood(){
		double logLikelihood = Math.log(sumScaledLikelihood) + logScaler;
		return logLikelihood;
	}

	public void scaleLogProb(double logProb){
		sumScaledLikelihood += scale(logProb, logScaler);
	}

	public void scaleLogProb(double logProb, int count){

		double sum = scale(logProb, logScaler) * count;
		sumScaledLikelihood += sum;
	}

	public void addScaledLogProb(double scaledLogProb){
		sumScaledLikelihood += scaledLogProb;
	}
	public double scale(double logProb){
		return scale(logProb, logScaler);
	}

	private static double scale(double logProb, double logScaler) {

		double expB = Math.exp(logProb - logScaler);
		return expB;
	}
	@Deprecated
	private ArrayList<Double> likelihood;
	@Deprecated
	public void addToList(double logProb){
		likelihood.add(scale(logProb, logScaler));
	}
	@Deprecated
	public void sumList(){
		double sumLikelihood = 0;
		for (Double d : likelihood) {
			sumLikelihood += d;
		}
		double logLikelihood = Math.log(sumLikelihood) + logScaler;
	}
}
