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
import srp.spectrum.operator.DeltaExchangeColumnSpectrumOperator;



public class DeltaExchangeColumnSpectrumOperatorTest {

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
		int spectrumCount = 4;
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, spectrumCount);
		DeltaExchangeColumnSpectrumOperator op = new DeltaExchangeColumnSpectrumOperator(
				spectrumModel, 0.1, CoercionMode.COERCION_OFF);

		double[][][] storedFrequencies = new double[spectrumCount][spectrumModel.getSiteCount()][4];
		for (int i = 0; i < storedFrequencies.length; i++) {
			Spectrum spectrum = spectrumModel.getSpectrum(i);
			for (int j = 0; j < storedFrequencies[i].length; j++) {
				storedFrequencies[i][j] = spectrum.getFrequenciesAt(j);
			}
		}

		for (int o = 0; o < 100; o++) {
			try {
				op.doOperation();
				
				SpectrumOperationRecord opRecord = spectrumModel.getSpectrumOperationRecord();

				int siteIndex = opRecord.getColumnIndex();
				double[] delta = opRecord.getDelta();
				
				for (int i = 0; i < spectrumCount; i++) {
					Spectrum spectrum = spectrumModel.getSpectrum(i);
					double[] frequencies = spectrum.getFrequenciesAt(siteIndex);
					int count = 0;
					for (int f = 0; f < frequencies.length; f++) {
						if(frequencies[f]!= storedFrequencies[i][siteIndex][f]){
							count++;
							double absDelta = Math.abs(frequencies[f]-storedFrequencies[i][siteIndex][f]);
							assertEquals(delta[i], absDelta, 1e-8);
						}
						storedFrequencies[i][siteIndex][f] = frequencies[f];
					}
					assertEquals(2, count);
				}
			
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}

	}

	

	@Test
	public void testDoOperator2() throws Exception {
		String[] seqs = new String[]{
				"AAACGTTT",
				"AAACGT..",
				"..AGGTTC",
				};
		int spectrumCount = 10;
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, spectrumCount);
		DeltaExchangeColumnSpectrumOperator op = new DeltaExchangeColumnSpectrumOperator(
				spectrumModel, 0.1, CoercionMode.COERCION_OFF);

		double[][][] storedFrequencies = new double[spectrumCount][spectrumModel.getSiteCount()][4];
		for (int i = 0; i < storedFrequencies.length; i++) {
			Spectrum spectrum = spectrumModel.getSpectrum(i);
			for (int j = 0; j < storedFrequencies[i].length; j++) {
				storedFrequencies[i][j] = spectrum.getFrequenciesAt(j);
			}
		}

		for (int o = 0; o < 10000; o++) {
			try {
				op.doOperation();
				
				SpectrumOperationRecord opRecord = spectrumModel.getSpectrumOperationRecord();

				int siteIndex = opRecord.getColumnIndex();
				double[] delta = opRecord.getDelta();
				
				for (int i = 0; i < spectrumCount; i++) {
					Spectrum spectrum = spectrumModel.getSpectrum(i);
					double[] frequencies = spectrum.getFrequenciesAt(siteIndex);
					int count = 0;
					for (int f = 0; f < frequencies.length; f++) {
						if(frequencies[f]!= storedFrequencies[i][siteIndex][f]){
							count++;
							double absDelta = Math.abs(frequencies[f]-storedFrequencies[i][siteIndex][f]);
							assertEquals(delta[i], absDelta, 1e-8);
						}
						storedFrequencies[i][siteIndex][f] = frequencies[f];
					}
					assertEquals(2, count);
				}
			
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}		
		
	}

}
