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
import srp.spectrum.operator.DirichletAlphaSpectrumOperator;
import srp.spectrum.operator.DirichletSpectrumOperator;
import test.TestUtils;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

public class DirichletAlphaSpectrumOperatorTest {

	public static final int DIMENSION = 4;
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
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 1, 2);
		DirichletAlphaSpectrumOperator op = new DirichletAlphaSpectrumOperator(
				spectrumModel, 100, CoercionMode.COERCION_OFF);

		double[][] storedFrequencies = new double[spectrumModel.getSiteCount()][DIMENSION];
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

				
				double[] frequencies = spectrumModel.getSpecturmFrequencies(spectrumIndex, siteIndex);
				
				int count = 0;
				double delta = 0;
				for (int f = 0; f < frequencies.length; f++) {
					if(frequencies[f]!= storedFrequencies[siteIndex][f]){
						count++;
						delta += (frequencies[f]-storedFrequencies[siteIndex][f]);
					}
					storedFrequencies[siteIndex][f] = frequencies[f];
				}
				assertEquals(0, delta, TestUtils.UNITTEST_THRESHOLD);
				assertEquals(4, count);
			
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}

	}
	


}
