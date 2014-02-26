package test.srp.spectrum.operator;


import static org.junit.Assert.assertEquals;

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
import srp.spectrum.operator.DeltaExchangeSingleSpectrumOperator;
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
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		int spectrumLength = aMap.getLength();
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
				
				SpectrumOperationRecord opRecord = spectrumModel.getSpectrumOperationRecord();
				int spectrumIndex = opRecord.getSpectrumIndex();
				int siteIndex = opRecord.getAllSiteIndexs()[0];
				double delta = opRecord.getDelta()[0];
				
				double[] frequencies = spectrumModel.getSpecturmFrequencies(spectrumIndex, siteIndex);
				
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
	public void testDoOperatorMultiSpectrum() throws Exception {
		String[] seqs = new String[]{
				"AAACGTTT",
				"AAACGT..",
				"..AGGTTC",
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		int spectrumLength = aMap.getLength();
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
				
				SpectrumOperationRecord opRecord = spectrumModel.getSpectrumOperationRecord();
				int spectrumIndex = opRecord.getSpectrumIndex();
				int siteIndex = opRecord.getAllSiteIndexs()[0];
				double delta = opRecord.getDelta()[0];
				
				Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
				double[] frequencies = spectrum.getFrequenciesAt(siteIndex);
				
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
