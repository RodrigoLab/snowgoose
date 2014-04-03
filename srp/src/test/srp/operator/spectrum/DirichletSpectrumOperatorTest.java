package test.srp.operator.spectrum;


import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.stat.StatUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.haplotypes.AlignmentUtils;
import srp.operator.spectrum.DirichletSpectrumOperator;
import srp.shortreads.AlignmentMapping;
import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperationRecord;
import test.TestUtils;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.GammaFunction;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;



public class DirichletSpectrumOperatorTest {

	public static final int DIMENSION = 4;
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}



	private double MIN_FREQ = SpectraParameter.MIN_FREQ;
	private double MAX_FREQ = SpectraParameter.MAX_FREQ;

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void testDoOperatorOneSpectrum() throws Exception {
		
		String[] seqs = new String[]{
				"AAACGTTT",
				"AAACGT..",
				"..AGGTTC",
				};
		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, 1);
		DirichletSpectrumOperator op = new DirichletSpectrumOperator(
				spectrumModel, 1, 100, CoercionMode.COERCION_OFF);

		double[][] storedFrequencies = new double[spectrumModel.getSpectrumLength()][DIMENSION];
		Spectrum spectrum = spectrumModel.getSpectrum(0);
		for (int i = 0; i < storedFrequencies.length; i++) {
			storedFrequencies[i] = spectrum.getFrequenciesAt(i);
		}
		
		for (int o = 0; o < 10000; o++) {
			try {
				op.doOperation();
				
				SpectrumOperationRecord opRecord = spectrumModel.getSpectrumOperationRecord();
				int spectrumIndex = opRecord.getSpectrumIndex();
				int siteIndex = opRecord.getAllSiteIndexs()[0];
//				double[] delta = opRecord.getDelta();
				
				double[] frequencies = spectrumModel.getSpecturmFrequencies(spectrumIndex, siteIndex);
				
				int count = 0;
				double delta = 0;
				double absDelta = 0;
//					System.out.println(Arrays.toString(frequencies));
//					System.out.println(Arrays.toString(storedFrequencies[s]));
				for (int f = 0; f < frequencies.length; f++) {
					if(frequencies[f]!= storedFrequencies[siteIndex][f]){
						count++;
						absDelta += Math.abs(frequencies[f]-storedFrequencies[siteIndex][f]);
						delta += (frequencies[f]-storedFrequencies[siteIndex][f]);
//							assertEquals(delta[i], absDelta, 1e-8);
					}
					else if(frequencies[f]==MIN_FREQ | frequencies[f]==MAX_FREQ){
						count++;
					}
					storedFrequencies[siteIndex][f] = frequencies[f];
				}
//					System.out.println(delta +"\t"+ absDelta);
//					System.out.println();
				assertEquals(0, delta, 1-MAX_FREQ);
				assertEquals(4, count);
			
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}

	}
	

	@Test
	public void testProposalRatio() throws Exception {
		
		String[] seqs = new String[]{
				"AAACGTTT",
				"AAACGT..",
				"..AGGTTC",
				};
		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, 5);
		DirichletSpectrumOperator op = new DirichletSpectrumOperator(
				spectrumModel, 1, 100, CoercionMode.COERCION_OFF);
		double alpha = op.getAlpha();
		
		double[][][] storedFrequencies = new double[spectrumModel
				.getSpectrumCount()][spectrumModel.getSpectrumLength()][DIMENSION];

		for (int s = 0; s < storedFrequencies.length; s++) {
			Spectrum spectrum = spectrumModel.getSpectrum(s);
			for (int l = 0; l < storedFrequencies[s].length; l++) {
				storedFrequencies[s][l] = spectrum.getFrequenciesAt(l);
			}
		}
		for (int o = 0; o < 10000; o++) {
			try {
				double ratio = op.doOperation();
				
				SpectrumOperationRecord opRecord = spectrumModel.getSpectrumOperationRecord();
				int spectrumIndex = opRecord.getSpectrumIndex();
				int siteIndex = opRecord.getAllSiteIndexs()[0];
				
				SpectraParameter spectra = spectrumModel.getSpectrum(spectrumIndex).getSpectra(siteIndex);
				
				double[] newFreq = spectra.getFrequencies();
				double[] oldFreq = storedFrequencies[spectrumIndex][siteIndex];
				double expected = calculateLogq(alpha, oldFreq, newFreq);
				assertEquals(expected, ratio, TestUtils.UNITTEST_THRESHOLD);
				System.arraycopy(newFreq, 0, storedFrequencies[spectrumIndex][siteIndex], 0, DIMENSION);
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}
	}

	private static double calculateLogq(double alphaPi, double[] oldPi, double[] newPi){
		
		int nStates = DIMENSION;
		int i;
		double sum = 0.0;
		double x, y;
		//Copied from MrBayes mcmc.c line 34134, same as line 33945
		sum = 0.0;
		for (i=0; i<nStates; i++)
			sum += newPi[i]*alphaPi;
		x = LnGamma(sum);
		for (i=0; i<nStates; i++)
			x -= LnGamma(newPi[i]*alphaPi);
		for (i=0; i<nStates; i++)
			x += (newPi[i]*alphaPi-1.0)*log(oldPi[i]);
		sum = 0.0;
		for (i=0; i<nStates; i++)
			sum += oldPi[i]*alphaPi;
		y = LnGamma(sum);
		for (i=0; i<nStates; i++)
			y -= LnGamma(oldPi[i]*alphaPi);
		for (i=0; i<nStates; i++)
			y += (oldPi[i]*alphaPi-1.0)*log(newPi[i]);
		//
		double ratio = x - y;

		return ratio;
	}
	private static double LnGamma(double x){
		return GammaFunction.lnGamma(x);
	}
	private static double log(double x){
		return Math.log(x);
	}
			
/*
	 get proposal ratio 
	sum = 0.0;
	for (i=0; i<nStates; i++)
		sum += newPi[i]*alphaPi;
	x = LnGamma(sum);
	for (i=0; i<nStates; i++)
		x -= LnGamma(newPi[i]*alphaPi);
	for (i=0; i<nStates; i++)
		x += (newPi[i]*alphaPi-1.0)*log(oldPi[i]);
	sum = 0.0;
	for (i=0; i<nStates; i++)
		sum += oldPi[i]*alphaPi;
	y = LnGamma(sum);
	for (i=0; i<nStates; i++)
		y -= LnGamma(oldPi[i]*alphaPi);
	for (i=0; i<nStates; i++)
		y += (oldPi[i]*alphaPi-1.0)*log(newPi[i]);
	(*lnProposalRatio) = x - y;
	
	
	//mb.h:#define PI_MIN				    0.000001f

	    do {
        DirichletRandomVariable (dirichletParameters, newPi, nStates, seed);
        isValid = YES;
        for (i=0; i<nStates; i++)
            {
		    if (newPi[i] < PI_MIN)
                {
                isValid = NO;
                break;
                }
            }
		} while (isValid == NO);
		
///////////////
 
 Autotune Dirichlet move 
void AutotuneDirichlet (MrBFlt acceptanceRate, MrBFlt targetRate, int batch, MrBFlt *alphaPi, MrBFlt minTuning, MrBFlt maxTuning)
{
    MrBFlt delta, logTuning, newTuning;

    delta = 1.0 / sqrt(batch);
    delta = 0.01 < delta ? 0.01 : delta;

    logTuning = log(*alphaPi);

    if (acceptanceRate > targetRate)
        logTuning -= delta;
    else
        logTuning += delta;

    newTuning = exp(logTuning);
    if (newTuning > minTuning && newTuning < maxTuning)
        *alphaPi = newTuning;
}

		
	*/
	

	
	@Test
	public void testDoOperatorMultiSpectrum() throws Exception {
		String[] seqs = new String[]{
				"AAACGTTT",
				"AAACGT..",
				"..AGGTTC",
				};
		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, 5);
		DirichletSpectrumOperator op = new DirichletSpectrumOperator(
				spectrumModel,  10, 100, CoercionMode.COERCION_OFF);

		double[][][] storedFrequencies = new double[spectrumModel
				.getSpectrumCount()][spectrumModel.getSpectrumLength()][4];
		for (int s = 0; s < storedFrequencies.length; s++) {
			
			Spectrum spectrum = spectrumModel.getSpectrum(s);
			for (int l = 0; l < storedFrequencies[s].length; l++) {
				storedFrequencies[s][l] = spectrum.getFrequenciesAt(l);
			}
		}
		for (int o = 0; o < 10000; o++) {
			try {
				op.doOperation();
				
				SpectrumOperationRecord opRecord = spectrumModel.getSpectrumOperationRecord();
				int spectrumIndex = opRecord.getSpectrumIndex();
				int[] siteIndexs = opRecord.getAllSiteIndexs();
//				double[] delta = opRecord.getDelta();
				
				Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
				for (int i = 0; i < siteIndexs.length; i++) {
					
					double[] frequencies = spectrum.getFrequenciesAt(siteIndexs[i]);
					int count = 0;
					double delta = 0;
					double absDelta = 0;
					double[] spectraFrequencies = storedFrequencies[spectrumIndex][siteIndexs[i]];
					for (int f = 0; f < frequencies.length; f++) {
						if(frequencies[f]!= spectraFrequencies[f]){
							count++;
							delta += (frequencies[f]-spectraFrequencies[f]);
							absDelta += Math.abs(frequencies[f]-spectraFrequencies[f]);
//							assertEquals(delta[i], absDelta, 1e-8);
						}
						else if(frequencies[f]==MIN_FREQ | frequencies[f]==MAX_FREQ){
							count++;
						}
						spectraFrequencies[f] = frequencies[f];
					}
					assertEquals(0, delta, 1-MAX_FREQ);
					assertEquals(4, count);
				}
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}		
		
	}
//RevBayes - simplexMove.cpp & DistributionDirichlet.cpp
//MrBayes code
	/*
		 get proposal ratio 
		sum = 0.0;
		for (i=0; i<nStates; i++)
			sum += newPi[i]*alphaPi;
		x = LnGamma(sum);
		for (i=0; i<nStates; i++)
			x -= LnGamma(newPi[i]*alphaPi);
		for (i=0; i<nStates; i++)
			x += (newPi[i]*alphaPi-1.0)*log(oldPi[i]);
		sum = 0.0;
		for (i=0; i<nStates; i++)
			sum += oldPi[i]*alphaPi;
		y = LnGamma(sum);
		for (i=0; i<nStates; i++)
			y -= LnGamma(oldPi[i]*alphaPi);
		for (i=0; i<nStates; i++)
			y += (oldPi[i]*alphaPi-1.0)*log(newPi[i]);
		(*lnProposalRatio) = x - y;
		
		
		//mb.h:#define PI_MIN				    0.000001f
	
		    do {
	        DirichletRandomVariable (dirichletParameters, newPi, nStates, seed);
	        isValid = YES;
	        for (i=0; i<nStates; i++)
	            {
			    if (newPi[i] < PI_MIN)
	                {
	                isValid = NO;
	                break;
	                }
	            }
			} while (isValid == NO);
			
	///////////////
	 
	 Autotune Dirichlet move 
	void AutotuneDirichlet (MrBFlt acceptanceRate, MrBFlt targetRate, int batch, MrBFlt *alphaPi, MrBFlt minTuning, MrBFlt maxTuning)
	{
	    MrBFlt delta, logTuning, newTuning;
	
	    delta = 1.0 / sqrt(batch);
	    delta = 0.01 < delta ? 0.01 : delta;
	
	    logTuning = log(*alphaPi);
	
	    if (acceptanceRate > targetRate)
	        logTuning -= delta;
	    else
	        logTuning += delta;
	
	    newTuning = exp(logTuning);
	    if (newTuning > minTuning && newTuning < maxTuning)
	        *alphaPi = newTuning;
	}
	
			
		*/
		
	
		
		public void testRandomGammaSpeed() throws Exception {
		
			double alpha = 0.5;
	//		DirichletDistribution dd = new DirichletDistribution(alphas);
			GammaDistribution drGamma = new GammaDistribution(alpha, 1);
			long time1 = System.currentTimeMillis();
			for (int t = 0; t < 1e6; t++) {alpha = MathUtils.nextDouble()*2;
				drGamma.nextGamma();
			}
			long time2 = System.currentTimeMillis();
			System.out.println("Single drGamma:\t"+(time2 - time1) + "\t");
			alpha = MathUtils.nextDouble()*2;
			time1 = System.currentTimeMillis();
			for (int t = 0; t < 1e6; t++) {alpha = MathUtils.nextDouble()*2;
				drGamma.setScale(alpha);
				drGamma.nextGamma();
			}
			time2 = System.currentTimeMillis();
			System.out.println("Multi drGamma:\t"+(time2 - time1) + "\t");
			
			
			time1 = System.currentTimeMillis();
			for (int t = 0; t < 1e6; t++) {alpha = MathUtils.nextDouble()*2;
				alpha= MathUtils.nextGamma(alpha, 1);
			}
			System.out.println(alpha);
			time2 = System.currentTimeMillis();
			System.out.println("Multi MuthUtil:\t"+(time2 - time1) + "\t");
			
	
			RandomDataGenerator rd = new RandomDataGenerator();
			time1 = System.currentTimeMillis();
			for (int t = 0; t < 1e6; t++) {alpha = MathUtils.nextDouble()*2;
				rd.nextGamma(alpha, 1);
			}
			time2 = System.currentTimeMillis();
			System.out.println("Multi rdg:\t"+(time2 - time1) + "\t");
			
	//		org.apache.commons.math3.distribution.GammaDistribution comGamma = new org.apache.commons.math3.distribution.GammaDistribution(alpha, 1);
	//		time1 = System.currentTimeMillis();
	//		for (int t = 0; t < 1e6; t++) {
	//			comGamma = new org.apache.commons.math3.distribution.GammaDistribution(alpha, 1);
	//			comGamma.sample();
	//		}
	//		time2 = System.currentTimeMillis();
	//		System.out.println("Multi drGamma:\t"+(time2 - time1) + "\t");
			
	//		alpha = 0.5
	//		Single drGamma:	1961	
	//		Multi drGamma:	1261	
	//		Multi MuthUtil:	487
	//		Multi rdg:	2287	
	//		Multi drGamma:	28939
	//		
	//		alpha = 1.5
	//		Single drGamma:	2236	
	//		Multi drGamma:	2198	
	//		Multi MuthUtil:	380	
	//		Multi rdg:	1631	
	//		Multi drGamma:	32053	
	
	
		}

}
