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
import srp.operator.spectrum.DeltaExchangeColumnSpectrumOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;



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
	public void testConstructor() throws Exception {
		String[] seqs = new String[]{
				"AAAC",
				"AACT",
				"ACGT"
				};
		int spectrumCount = 4;
		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, spectrumCount);
		DeltaExchangeColumnSpectrumOperator op = new DeltaExchangeColumnSpectrumOperator(
				spectrumModel, 0.1, CoercionMode.COERCION_OFF);
	}
	@Test
	public void testDoOperator() throws Exception {
		String[] seqs = new String[]{
				"AAAC",
				"AACT",
				"ACGT"
				};
		int spectrumCount = 4;
		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, spectrumCount);
		DeltaExchangeColumnSpectrumOperator op = new DeltaExchangeColumnSpectrumOperator(
				spectrumModel, 0.1, CoercionMode.COERCION_OFF);

		double[][][] storedFrequencies = new double[spectrumCount][spectrumModel.getSpectrumLength()][4];
		for (int i = 0; i < storedFrequencies.length; i++) {
			Spectrum spectrum = spectrumModel.getSpectrum(i);
			for (int j = 0; j < storedFrequencies[i].length; j++) {
				storedFrequencies[i][j] = spectrum.getFrequenciesAt(j);
			}
		}

		for (int o = 0; o < 100; o++) {
			try {
				op.doOperation();
				
				OperationRecord opRecord = spectrumModel.getOperationRecord();

				int siteIndex = opRecord.getSingleIndex();
//				double[] delta = opRecord.getDelta();
				
				for (int i = 0; i < spectrumCount; i++) {
					Spectrum spectrum = spectrumModel.getSpectrum(i);
					double[] frequencies = spectrum.getFrequenciesAt(siteIndex);
					int count = 0;
					double delta = 0;
					for (int f = 0; f < frequencies.length; f++) {
						if(frequencies[f]!= storedFrequencies[i][siteIndex][f]){
							count++;
							delta += (frequencies[f]-storedFrequencies[i][siteIndex][f]);
//							delta += Math.abs(frequencies[f]-storedFrequencies[i][siteIndex][f]);
						}

						storedFrequencies[i][siteIndex][f] = frequencies[f];
					}
					assertEquals(2, count);
					assertEquals(0, delta, 1e-8);
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
		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, spectrumCount);
		DeltaExchangeColumnSpectrumOperator op = new DeltaExchangeColumnSpectrumOperator(
				spectrumModel, 0.1, CoercionMode.COERCION_OFF);

		double[][][] storedFrequencies = new double[spectrumCount][spectrumModel.getSpectrumLength()][4];
		for (int i = 0; i < storedFrequencies.length; i++) {
			Spectrum spectrum = spectrumModel.getSpectrum(i);
			for (int j = 0; j < storedFrequencies[i].length; j++) {
				storedFrequencies[i][j] = spectrum.getFrequenciesAt(j);
			}
		}

		for (int o = 0; o < 10000; o++) {
			try {
				op.doOperation();
				
				OperationRecord opRecord = spectrumModel.getOperationRecord();

				int siteIndex = opRecord.getSingleIndex();
				
				
				for (int i = 0; i < spectrumCount; i++) {
					Spectrum spectrum = spectrumModel.getSpectrum(i);
					double[] frequencies = spectrum.getFrequenciesAt(siteIndex);
					int count = 0;
					double delta = 0;
					for (int f = 0; f < frequencies.length; f++) {
						if(frequencies[f]!= storedFrequencies[i][siteIndex][f]){
							count++;
							delta += (frequencies[f]-storedFrequencies[i][siteIndex][f]);
							
						}
						storedFrequencies[i][siteIndex][f] = frequencies[f];
					}
					assertEquals(2, count);
					assertEquals(0, delta, 1e-8);
				}
			
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}		
		
	}

}
