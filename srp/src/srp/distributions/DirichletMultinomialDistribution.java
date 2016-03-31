package srp.distributions;

import java.util.Arrays;

public class DirichletMultinomialDistribution {

	static final double PHI = 0.1;
	static final double ERROR_PROB = 0.0107;
	double alphas_total;
	
	double[][] alphas = new double[4][4];
	
	public DirichletMultinomialDistribution() {

		alphas_total = (1.0 - PHI) / PHI;
		double offDiagAlpha = ERROR_PROB / 3.0 * alphas_total;
		double diagAlpha = (1.0 - ERROR_PROB) * alphas_total;
		
		for (int i = 0; i < alphas.length; i++) {
			Arrays.fill(alphas[i], offDiagAlpha);
			alphas[i][i] = diagAlpha;
//			System.out.println(Arrays.toString(alphas[i]) );
		}

	}
	public double[] FourDirichletMultinomialLogProbability(int[] readCount){
		int readTotal = readCount[0] + readCount[1] + readCount[2] + readCount[3];
		double[] result = new double[4];
		for (int i = 0; i < result.length; i++) {
			result[i] = DirichletMultinomialLogProbability(i, readCount);
		}
		
		return result;
		
	}
	
	public double DirichletMultinomialLogProbability(int acgt,
			int[] readCount) {

		int readCountTotal = readCount[0] + readCount[1] + readCount[2] + readCount[3];

		
		double result = 0.0;
		for (int i = 0; i < readCount.length; i++) {
			for (int x = 0; x < readCount[i]; ++x) {
				result += Math.log(alphas[acgt][i] + x);
			}
		}
		for (int x = 0; x < readCountTotal; ++x) {
			result -= Math.log(alphas_total + x);
		}

		return result;
	}

}
