package test.srp.operator.spectrum;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.evolution.OperationRecord;
import srp.evolution.shortreads.AlignmentMapping;
import srp.evolution.spectrum.SpectraParameter;
import srp.evolution.spectrum.Spectrum;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import srp.haplotypes.AlignmentUtils;
import srp.operator.spectrum.AbstractDirichletSpectrumOperator;
import srp.operator.spectrum.DeltaExchangeSingleSpectrumOperator;
import srp.operator.spectrum.DirichletAlphaSpectrumOperator;
import srp.operator.spectrum.DirichletSpectrumOperator;
import test.TestUtils;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class DirichletAlphaSpectrumOperatorTest {

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
		AbstractDirichletSpectrumOperator op = new DirichletAlphaSpectrumOperator(
				spectrumModel, 100, CoercionMode.COERCION_OFF);

		double[][] storedFrequencies = new double[spectrumModel.getSpectrumLength()][DIMENSION];
		Spectrum spectrum = spectrumModel.getSpectrum(0);
		for (int i = 0; i < storedFrequencies.length; i++) {
			storedFrequencies[i] = spectrum.getFrequenciesAt(i);
		}
		
		for (int o = 0; o < 10000; o++) {
			try {
				int alpha = MathUtils.nextInt(1000)+1;
				op = new DirichletAlphaSpectrumOperator(
						spectrumModel, alpha, CoercionMode.COERCION_OFF);

				op.doOperation();
				
				OperationRecord opRecord = spectrumModel.getOperationRecord();
				int spectrumIndex = opRecord.getSpectrumIndex();
				int siteIndex = opRecord.getSingleIndex();

				
				double[] frequencies = spectrumModel.getSpecturmFrequencies(spectrumIndex, siteIndex);
				
				int count = 0;
				double delta = 0;
				for (int f = 0; f < frequencies.length; f++) {
					if(frequencies[f]!= storedFrequencies[siteIndex][f]){
						count++;
						delta += (frequencies[f]-storedFrequencies[siteIndex][f]);
					}
					else if(frequencies[f]==MIN_FREQ | frequencies[f]==MAX_FREQ){
						count++;
					}
					storedFrequencies[siteIndex][f] = frequencies[f];
				}
				assertEquals(Arrays.toString(storedFrequencies[siteIndex]) +"\n"+Arrays.toString(frequencies),
						0, delta, 1-MAX_FREQ);
				assertEquals(4, count);
			
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}

	}
	


}
