package test.srp.spectrum.operator;


import static org.junit.Assert.*;

import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperationRecord;
import srp.spectrum.likelihood.ShortReadsSpectrumLikelihood;
import srp.spectrum.operator.SingleSpectrumDeltaExchangeOperator;

public class SingleSpectrumDeltaExchangeOperatorTest {

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
	public void testDoOperator() throws Exception {
		String[] seqs = new String[]{
				"AAAC",
				"AACT",
				"ACGT"
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 1);
		SingleSpectrumDeltaExchangeOperator op = new SingleSpectrumDeltaExchangeOperator(
				spectrumModel, 0.1, CoercionMode.COERCION_OFF);

		double[][] storedFrequencies = new double[spectrumModel.getSiteCount()][4];
		Spectrum spectrum = spectrumModel.getSpectrum(0);
		for (int i = 0; i < storedFrequencies.length; i++) {
			storedFrequencies[i] = spectrum.getFrequencies(i);
		}
		
		for (int o = 0; o < 100; o++) {
			try {
				op.doOperation();
				
				SpectrumOperationRecord opRecord = spectrumModel.getSpectrumOperationRecord();
				int spectrumIndex = opRecord.getSpectrumIndex();
				int siteIndex = opRecord.getSiteIndex();
				double delta = opRecord.getDelta()[0];
				
				spectrum = spectrumModel.getSpectrum(spectrumIndex);
				double[] frequencies = spectrum.getFrequencies(siteIndex);
				
				int count = 0;
				for (int f = 0; f < frequencies.length; f++) {
					if(frequencies[f]!= storedFrequencies[siteIndex][f]){
						count++;
						double absDelta = Math.abs(frequencies[f]-storedFrequencies[siteIndex][f]);
						assertEquals(delta, absDelta, 1e-8);
					}
					storedFrequencies[siteIndex][f] = frequencies[f];
				}
				assertEquals(2, count);
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}

	}

	

	@Test
	public void testDoOperator2() throws Exception {
		String[] seqs = new String[]{
				"AAACGT",
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 5);
		SingleSpectrumDeltaExchangeOperator op = new SingleSpectrumDeltaExchangeOperator(
				spectrumModel, 0.1, CoercionMode.COERCION_OFF);

		double[][][] storedFrequencies = new double[spectrumModel
				.getSpectrumCount()][spectrumModel.getSiteCount()][4];
		for (int s = 0; s < storedFrequencies.length; s++) {
			
			Spectrum spectrum = spectrumModel.getSpectrum(s);
			for (int l = 0; l < storedFrequencies[s].length; l++) {
				storedFrequencies[s][l] = spectrum.getFrequencies(l);
			}
		}
		for (int o = 0; o < 10000; o++) {
			try {
				op.doOperation();
				
				SpectrumOperationRecord opRecord = spectrumModel.getSpectrumOperationRecord();
				int spectrumIndex = opRecord.getSpectrumIndex();
				int siteIndex = opRecord.getSiteIndex();
				double delta = opRecord.getDelta()[0];
				
				Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
				double[] frequencies = spectrum.getFrequencies(siteIndex);
				
				int count = 0;
				double[] spectraFrequencies = storedFrequencies[spectrumIndex][siteIndex];
				for (int f = 0; f < frequencies.length; f++) {
					if(frequencies[f]!= spectraFrequencies[f]){
						count++;
						double absDelta = Math.abs(frequencies[f]-spectraFrequencies[f]);
						assertEquals(delta, absDelta, 1e-8);
					}
					spectraFrequencies[f] = frequencies[f];
				}
				assertEquals(2, count);
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}		
		
	}

}
