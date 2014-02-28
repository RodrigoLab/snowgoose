package srp.haplotypes.likelihood;

import java.util.ArrayList;

import org.apache.commons.math3.util.FastMath;

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

	public void addLogProb(double logProb){
		sumScaledLikelihood += scale(logProb, logScaler);
	}

	public void minusLogProb(double logProb){
		sumScaledLikelihood -= scale(logProb, logScaler);
	}

	public void addLogProbMulti(double logProb, int count){
		double sum = scale(logProb, logScaler) * count;
		sumScaledLikelihood += sum;
	}


	public double scale(double logProb){
		return scale(logProb, logScaler);
//		return Math.exp(logProb - logScaler);
//		return (logProb - logScaler);
//		double a = Math.exp(logProb - logScaler);
//		double b = FastMath.exp(logProb - logScaler);
//		if(  ( Math.log(a) + logScaler  ) !=( Math.log(b) + logScaler  )  ){
//			System.out.println(a +"\t"+ b +"\t"+(a-b) +"\t"+ (Math.log( (a-b) ) + logScaler) +"\t"+  logProb);
//			System.out.println( (( Math.log(a) + logScaler  ) ==( Math.log(b) + logScaler  ) ) +"\t"+  ( Math.log(a) + logScaler  ) +"\t"+ ( Math.log(b) + logScaler  ) );
//			System.out.println();
//		}
//		return FastMath.exp(logProb - logScaler);
	}

	private static double scale(double logProb, double logScaler) {
//		double expB = Math.exp(logProb - logScaler);
		double expB = FastMath.exp(logProb - logScaler);

		return expB;
	}


	@Deprecated
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
