package test.srp.spectrum.operator;


import static org.junit.Assert.*;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;
import dr.math.distributions.DirichletDistribution;
import dr.math.distributions.GammaDistribution;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperationRecord;
import srp.spectrum.likelihood.ShortReadsSpectrumLikelihood;
import srp.spectrum.operator.DirichletSpectrumOperator;



public class DirichletSpectrumOperatorTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

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
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 1);
		DirichletSpectrumOperator op = new DirichletSpectrumOperator(
				spectrumModel, 1, CoercionMode.COERCION_OFF);

//		double[][] storedFrequencies = new double[spectrumModel.getSiteCount()][4];
//		Spectrum spectrum = spectrumModel.getSpectrum(0);
//		for (int i = 0; i < storedFrequencies.length; i++) {
//			storedFrequencies[i] = spectrum.getFrequenciesAt(i);
//		}
		
		for (int o = 0; o < 100; o++) {
			try {
				op.doOperation();
				
//				SpectrumOperationRecord opRecord = spectrumModel.getSpectrumOperationRecord();
//				int spectrumIndex = opRecord.getSpectrumIndex();
//				int[] siteIndexs = opRecord.getAllSiteIndexs();
//				double[] delta = opRecord.getDelta();
//				
//				spectrum = spectrumModel.getSpectrum(spectrumIndex);
//				for (int i = 0; i < siteIndexs.length; i++) {
//					
//					double[] frequencies = spectrum.getFrequenciesAt(siteIndexs[i]);
//					int count = 0;
//					for (int f = 0; f < frequencies.length; f++) {
//						if(frequencies[f]!= storedFrequencies[siteIndexs[i]][f]){
//							count++;
//							double absDelta = Math.abs(frequencies[f]-storedFrequencies[siteIndexs[i]][f]);
//							assertEquals(delta[i], absDelta, 1e-8);
//						}
//						storedFrequencies[siteIndexs[i]][f] = frequencies[f];
//					}
//					assertEquals(2, count);
//				}
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}

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
	

	
	//	@Test
	public void testDoOperatorMultiSpectrum() throws Exception {
		String[] seqs = new String[]{
				"AAACGTTT",
				"AAACGT..",
				"..AGGTTC",
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 5);
		DirichletSpectrumOperator op = new DirichletSpectrumOperator(
				spectrumModel,  10, CoercionMode.COERCION_OFF);

		double[][][] storedFrequencies = new double[spectrumModel
				.getSpectrumCount()][spectrumModel.getSiteCount()][4];
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
				double[] delta = opRecord.getDelta();
				
				Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
				for (int i = 0; i < siteIndexs.length; i++) {
					
					double[] frequencies = spectrum.getFrequenciesAt(siteIndexs[i]);
					int count = 0;
					double[] spectraFrequencies = storedFrequencies[spectrumIndex][siteIndexs[i]];
					for (int f = 0; f < frequencies.length; f++) {
						if(frequencies[f]!= spectraFrequencies[f]){
							count++;
							double absDelta = Math.abs(frequencies[f]-spectraFrequencies[f]);
							assertEquals(delta[i], absDelta, 1e-8);
						}
						spectraFrequencies[f] = frequencies[f];
					}
					assertEquals(2, count);
				}
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}		
		
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
