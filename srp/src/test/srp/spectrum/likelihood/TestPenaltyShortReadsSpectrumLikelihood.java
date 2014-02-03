package test.srp.spectrum.likelihood;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.apache.commons.math3.stat.StatUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.likelihood.ShortReadsSpectrumLikelihood;
import dr.evolution.alignment.SimpleAlignment;

public class TestPenaltyShortReadsSpectrumLikelihood {

	public static final double ERROR = ShortReadsSpectrumLikelihood.ERROR_RATE;
	public static final double NOT_ERROR = ShortReadsSpectrumLikelihood.NOT_ERROR_RATE;
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
	public void testCalculateLikelihoodIdentical() {

		String[] seqs = new String[]{"AACCGGTT"};
	
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
//		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		
		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, alignment);
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double logLikelihood = likelihood.getLogLikelihood();
		double expected = Math.log(1*NOT_ERROR+0*ERROR)*8;
//		System.out.println((0.25*NOT_ERROR+0.75*ERROR) +"\t"+ Math.log(1*NOT_ERROR+0*ERROR) );
		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);
		
		
		String[] seqs2 = new String[]{"AACCGGTT","AACCGGTT"};
		alignment = AlignmentUtils.createAlignment(seqs2);
		spectrumModel = new SpectrumAlignmentModel(aMap, alignment);
		likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double logLikelihood2 = likelihood.getLogLikelihood();
		double expected2 = Math.log(1*NOT_ERROR+0*ERROR)*8+Math.log(2);
		assertEquals("0 mismatch",expected2, logLikelihood2, 1e-10);
		
		String[] seqs3 = new String[]{"AACCGGTT","AACCGGTT","AACCGGTT"};
		alignment = AlignmentUtils.createAlignment(seqs3);
		spectrumModel = new SpectrumAlignmentModel(aMap, alignment);
		likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double logLikelihood3 = likelihood.getLogLikelihood();
		double expected3 = Math.log(1*NOT_ERROR+0*ERROR)*8+Math.log(3);
		assertEquals("0 mismatch",expected3, logLikelihood3, 1e-10);
		
		double p = 0;
		p = Math.log(8);
//		p = Math.log(1*NOT_ERROR+0*ERROR)*8; 
//		p = Math.log(0.5*NOT_ERROR+0.5*ERROR)*8;
//		p = Math.log(0*NOT_ERROR+1*ERROR)*8;
		
		double expected1a = p+Math.log(1)*2;
		double expected2a = p+Math.log(2)*2;
		double expected3a = p+Math.log(3)*2;
//		
//		double expected1a = p+1*2;
//		double expected2a = p+2*2;
//		double expected3a = p+3*2;
//		
		System.out.println();
		System.out.println(expected +"\t"+ (expected-expected1a));
		System.out.println(expected2 +"\t"+ (expected2-expected2a));
		System.out.println(expected3 +"\t"+ (expected3-expected3a));
		
	}
	
	@Test
	public void testCalculateLikelihoodFixedSpectrum() {
		String[] seqs = new String[]{
				".AA",
				".AC",
				".GT"
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		seqs = new String[]{"AAA"};
		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, alignment);
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double loglLikelihood = likelihood.getLogLikelihood();
		double[] eachLikelihood = likelihood.unittestMethodGetEachLikelihood();
		
//		System.out.println(Arrays.toString(eachLikelihood));
		double[] expecteds = new double[]{ 
				0+Math.log(1*NOT_ERROR+0*ERROR)*2,
				0+Math.log(1*NOT_ERROR+0*ERROR)*1+Math.log(0*NOT_ERROR+1*ERROR)*1,
				0+Math.log(1*NOT_ERROR+0*ERROR)*0+Math.log(0*NOT_ERROR+1*ERROR)*2
			};
		assertArrayEquals(expecteds, eachLikelihood, 1e-8);
		assertEquals(StatUtils.sum(expecteds), loglLikelihood, 1e-8);
		
		
		
		
		
		seqs = new String[]{"AAA", "AAA"};
		alignment = AlignmentUtils.createAlignment(seqs);
		spectrumModel = new SpectrumAlignmentModel(aMap, alignment);
		likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double loglLikelihood2 = likelihood.getLogLikelihood();
		double[] eachLikelihood2 = likelihood.unittestMethodGetEachLikelihood();
//		System.out.println(Arrays.toString(eachLikelihood));
		double[] expecteds2 = new double[]{ 
				0+Math.log(1*NOT_ERROR+0*ERROR)*2+Math.log(2),
				0+Math.log(1*NOT_ERROR+0*ERROR)*1+Math.log(0*NOT_ERROR+1*ERROR)*1+Math.log(2),
				0+Math.log(1*NOT_ERROR+0*ERROR)*0+Math.log(0*NOT_ERROR+1*ERROR)*2+Math.log(2)
			};
		assertArrayEquals(expecteds2, eachLikelihood2, 1e-8);
		assertEquals(StatUtils.sum(expecteds2), loglLikelihood2, 1e-8);
		
		double p = 0;
		p = Math.log(3);
		p = Math.log(0.5*NOT_ERROR+0.5*ERROR)*3;
		
		double expected1a = p+Math.log(1)*2;
		double expected2a = p+Math.log(2)*2;

		System.out.println();
		System.out.println(loglLikelihood +"\t"+ (loglLikelihood-expected1a));
		System.out.println(loglLikelihood2 +"\t"+ (loglLikelihood2-expected2a));
		
		
	}

	
	@Test
	public void testCalculateLikelihoodSpectrum() {
		String[] seqs = new String[]{
				".AA",
				".AC",
				".GT"
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
			
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 1, 0);
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double loglLikelihood = likelihood.getLogLikelihood();
		double[] eachLikelihood = likelihood.unittestMethodGetEachLikelihood();
//		System.out.println(Arrays.toString(eachLikelihood));
		double[] expecteds = new double[]{ 
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*2,
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*1+Math.log(0.25*NOT_ERROR+0.75*ERROR)*1,
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*2
			};
		assertArrayEquals(expecteds, eachLikelihood, 1e-8);
		assertEquals(StatUtils.sum(expecteds), loglLikelihood, 1e-8);
		

		
		spectrumModel = new SpectrumAlignmentModel(aMap, 2, 0);
		likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double loglLikelihood2 = likelihood.getLogLikelihood();
		double[] eachLikelihood2 = likelihood.unittestMethodGetEachLikelihood();
//		System.out.println(Arrays.toString(eachLikelihood));
		double[] expecteds2 = new double[]{ 
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*2+Math.log(2),
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*1+Math.log(0.25*NOT_ERROR+0.75*ERROR)*1+Math.log(2),
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*2+Math.log(2)
			};
		assertArrayEquals(expecteds2, eachLikelihood2, 1e-8);
		assertEquals(StatUtils.sum(expecteds2), loglLikelihood2, 1e-8);
		

		double p = 0;
		p = Math.log(3);
		p = Math.log(0.5*NOT_ERROR+0.5*ERROR)*3;
		
		double expected1a = p+Math.log(1)*2;
		double expected2a = p+Math.log(2)*2;

		System.out.println();
		System.out.println(loglLikelihood +"\t"+ (loglLikelihood-expected1a));
		System.out.println(loglLikelihood2 +"\t"+ (loglLikelihood2-expected2a));
		
		
		
	}



	@Test
	public void testCalculateLikelihoodCustomSpectrum() {
		String[] seqs = new String[]{
				"AAAC",
				"AACT",
				"ACGT"
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 1, 1);
		Spectrum spectrum = spectrumModel.getSpectrum(0);
		for (int i = 0; i < spectrum.getLength(); i++) {
			double[] freqs = new double[]{1-(0.1*i*3), 0.1*i, 0.1*i, 0.1*i};
			spectrum.resetFrequencies(i, freqs);
//			System.out.println("SITE: "+i +"\t"+  Arrays.toString(spectrum.getFrequencies(i)));
		}
//		spectrumModel.setSpectrum(0, spectrum);
		
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double[] eachLikelihood = likelihood.unittestMethodGetEachLikelihood();
//		System.out.println(Arrays.toString(eachLikelihood));
		//Site1: 1, 0, 0, 0
		//Site2: 0.7, 0.1, 0.1, 0.1
		//Site3: 0.4, 0.2, 0.2, 0.2
		//Site4: 0.1, 0.3, 0.3, 0.3
		double[] expecteds = new double[]{ 
				//MMMD
				Math.log((1*NOT_ERROR+0*ERROR)*(0.7*NOT_ERROR+0.3*ERROR)*
						(0.4*NOT_ERROR+0.6*ERROR)*(0.3*NOT_ERROR+0.7*ERROR)),
				
				//MMDD
				Math.log(1*NOT_ERROR+0*ERROR)+Math.log(0.7*NOT_ERROR+0.3*ERROR)+
				Math.log(0.2*NOT_ERROR+0.8*ERROR)+Math.log(0.3*NOT_ERROR+0.7*ERROR),
				
				//MDDD
				Math.log(1*NOT_ERROR+0*ERROR)+Math.log(0.1*NOT_ERROR+0.9*ERROR)+
				Math.log(0.2*NOT_ERROR+0.8*ERROR)+Math.log(0.3*NOT_ERROR+0.7*ERROR),

			};
		assertArrayEquals(expecteds, eachLikelihood, 1e-8);
		

	}

}
