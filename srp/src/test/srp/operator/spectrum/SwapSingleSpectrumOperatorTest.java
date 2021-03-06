package test.srp.operator.spectrum;


import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.evolution.OperationRecord;
import srp.evolution.spectrum.Spectrum;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import srp.operator.spectrum.SwapSingleSpectrumOperator;
import dr.inference.operators.OperatorFailedException;

public class SwapSingleSpectrumOperatorTest {

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
		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, 5);
		SwapSingleSpectrumOperator op = new SwapSingleSpectrumOperator(spectrumModel);

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
				
				OperationRecord opRecord = spectrumModel.getOperationRecord();
				int spectrumIndex = opRecord.getSpectrumIndex();
				int siteIndex = opRecord.getSingleIndex();

				double[] frequencies = spectrumModel.getSpecturmFrequencies(spectrumIndex, siteIndex);
				double[] spectraFrequencies = storedFrequencies[spectrumIndex][siteIndex];

				int count = 0;
				int[] matches = new int[4];
				for (int f1 = 0; f1 < frequencies.length; f1++) {
					for (int f2 = 0; f2 < spectraFrequencies.length; f2++) {
						if(frequencies[f1]==spectraFrequencies[f2]){
							matches[f1]=f2;
							if(f1!=f2){
								count++;
							}
							break;
						}
					}
				}
				assertEquals(2, count);
				for (int f1 = 0; f1 < matches.length; f1++) {
					if(f1!=matches[f1]){
//						System.out.println(Arrays.toString(frequencies));
//						System.out.println(Arrays.toString(spectraFrequencies));
//						System.out.println(Arrays.toString(matches));
//						System.out.println();
						assertEquals(f1, matches[matches[f1]]);
					}
				}
				
				storedFrequencies[spectrumIndex][siteIndex] = frequencies;
				
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}		

		
		
	}

}
