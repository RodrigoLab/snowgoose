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
import srp.operator.spectrum.DeltaExchangeSingleSpectrumOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

public class DeltaExchangeSingleSpectrumOperatorTest {

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
	public void testDoOperation() throws Exception {
		String[] seqs = new String[]{
				"AAACGTTT",
				"AAACGT..",
				"..AGGTTC",
				};

		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, 1);
		DeltaExchangeSingleSpectrumOperator op = new DeltaExchangeSingleSpectrumOperator(
				spectrumModel, 0.1, CoercionMode.COERCION_OFF);

		double[][] storedFrequencies = new double[spectrumModel.getSpectrumLength()][4];
		Spectrum spectrum = spectrumModel.getSpectrum(0);
		for (int i = 0; i < storedFrequencies.length; i++) {
			storedFrequencies[i] = spectrum.getFrequenciesAt(i);
		}
		
		for (int o = 0; o < 100; o++) {
			try {
				op.doOperation();
				
				OperationRecord opRecord = spectrumModel.getOperationRecord();
				int spectrumIndex = opRecord.getSpectrumIndex();
				int siteIndex = opRecord.getSingleIndex();
//				double delta = opRecord.getDelta()[0];
				
				double[] frequencies = spectrumModel.getSpecturmFrequencies(spectrumIndex, siteIndex);
				
				int count = 0;
				double delta = 0;
				for (int f = 0; f < frequencies.length; f++) {
					if(frequencies[f]!= storedFrequencies[siteIndex][f]){
						count++;
						delta += (frequencies[f]-storedFrequencies[siteIndex][f]);
//						assertEquals(delta, absDelta, 1e-8);
					}
					storedFrequencies[siteIndex][f] = frequencies[f];
				}
				assertEquals(2, count);
				assertEquals(0, delta, 1e-8);
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}

	}

	

	@Test
	public void testDoOperatorMultiSpectrum() throws Exception {
		String[] seqs = new String[]{
				"AAACGTTT",
				"AAACGT..",
				"..AGGTTC",
				};
		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, 5);
		DeltaExchangeSingleSpectrumOperator op = new DeltaExchangeSingleSpectrumOperator(
				spectrumModel, 0.1, CoercionMode.COERCION_OFF);

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
//				double delta = opRecord.getDelta()[0];
				
				Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
				double[] frequencies = spectrum.getFrequenciesAt(siteIndex);
				
				int count = 0;
				double delta = 0;
				double[] spectraFrequencies = storedFrequencies[spectrumIndex][siteIndex];
				for (int f = 0; f < frequencies.length; f++) {
					if(frequencies[f]!= spectraFrequencies[f]){
						count++;
						delta += (frequencies[f]-spectraFrequencies[f]);
//						assertEquals(delta, absDelta, 1e-8);
					}
					spectraFrequencies[f] = frequencies[f];
				}
				assertEquals(2, count);
				assertEquals(0, delta, 1e-8);
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}		
		
	}


}
