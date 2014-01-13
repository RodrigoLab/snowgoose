package test.srp.spectrum.operator;


import static org.junit.Assert.*;
import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperationRecord;
import srp.spectrum.operator.MultiSpectrumDeltaExchangeOperator;
import srp.spectrum.operator.RecombinationSpectrumOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

public class RecombinationSpectrumOperatorTest {

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
	public void testOperator() throws Exception {
		String[] seqs = new String[]{
				"AAACGTTT",
				"AAACGT..",
				"..AGGTTC",
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 5);
		RecombinationSpectrumOperator op = new RecombinationSpectrumOperator(
				spectrumModel, 5, CoercionMode.COERCION_OFF);

		double[][][] storedFrequencies = new double[spectrumModel
				.getSpectrumCount()][spectrumModel.getSiteCount()][4];
		for (int s = 0; s < storedFrequencies.length; s++) {
			Spectrum spectrum = spectrumModel.getSpectrum(s);
			for (int l = 0; l < storedFrequencies[s].length; l++) {
				storedFrequencies[s][l] = spectrum.getFrequenciesAt(l);
			}
		}
		
		int spectrumLength = spectrumModel.getSpectrumLength();
		for (int o = 0; o < 10000; o++) {
			try {
				op.doOperation();
				
				SpectrumOperationRecord opRecord = spectrumModel.getSpectrumOperationRecord();
				
				int[] twoSpectrums = opRecord.getRecombinationSpectrumIndex();
				int[] twoPositions = opRecord.getRecombinationPositionIndex();
				System.out.println(Arrays.toString(twoPositions));
				Spectrum spectrum0 = spectrumModel.getSpectrum(twoSpectrums[0]);
				Spectrum spectrum1 = spectrumModel.getSpectrum(twoSpectrums[1]);
				
				for (int i = 0; i < spectrumLength; i++) {
					if(twoPositions[0] < i && i < twoPositions[1]){
						double[] frequencies = spectrum0.getFrequenciesAt(i);
						double[] expecteds = storedFrequencies[twoSpectrums[1]][i];
						assertArrayEquals(expecteds, frequencies, 0);
						
						frequencies = spectrum1.getFrequenciesAt(i);
						expecteds = storedFrequencies[twoSpectrums[0]][i];
						assertArrayEquals(expecteds, frequencies, 0);
					}
					else{
						double[] frequencies = spectrum0.getFrequenciesAt(i);
						double[] expecteds = storedFrequencies[twoSpectrums[0]][i];
						assertArrayEquals(expecteds, frequencies, 0);
						
						frequencies = spectrum1.getFrequenciesAt(i);
						expecteds = storedFrequencies[twoSpectrums[1]][i];
						assertArrayEquals(expecteds, frequencies, 0);
						
					}
					
				}
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}		

		
		
	}

}
