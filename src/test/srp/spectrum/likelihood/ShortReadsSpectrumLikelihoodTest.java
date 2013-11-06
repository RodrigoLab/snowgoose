package test.srp.spectrum.likelihood;


import static org.junit.Assert.*;

import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.HaplotypeModel;
import srp.likelihood.ShortReadLikelihood;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.likelihood.ShortReadsSpectrumLikelihood;
import dr.evolution.alignment.SimpleAlignment;

public class ShortReadsSpectrumLikelihoodTest {

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
		System.out.println((0.25*NOT_ERROR+0.75*ERROR) +"\t"+ Math.log(1*NOT_ERROR+0*ERROR) );
		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);
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
		
		double[] eachLikelihood = likelihood.getEachLikelihood();
		System.out.println(Arrays.toString(eachLikelihood));
		double[] expecteds = new double[]{ 
				0+Math.log(1*NOT_ERROR+0*ERROR)*2,
				0+Math.log(1*NOT_ERROR+0*ERROR)*1+Math.log(0*NOT_ERROR+1*ERROR)*1,
				0+Math.log(1*NOT_ERROR+0*ERROR)*0+Math.log(0*NOT_ERROR+1*ERROR)*2
			};
		assertArrayEquals(expecteds, eachLikelihood, 1e-8);
	}

	
	@Test
	public void testCalculateLikelihoodSpectrum() {
			String[] seqs = new String[]{
					".AA",
					".AC",
					".GT"
					};
			AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
			
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 1);
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double[] eachLikelihood = likelihood.getEachLikelihood();
		System.out.println(Arrays.toString(eachLikelihood));
		double[] expecteds = new double[]{ 
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*2,
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*1+Math.log(0.25*NOT_ERROR+0.75*ERROR)*1,
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*2
			};
		assertArrayEquals(expecteds, eachLikelihood, 1e-8);
		
//		logLikelihood = likelihood .getLogLikelihood();
//		expected = -0.086061253223681313806*4; //dbinom(0,8,E,log=T)
//		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);

//		double[] expecteds = new double[]{ 
//				0+Math.log((0.25*NOT_ERROR+0.75*ERROR)*2),
//				0+Math.log((0.25*NOT_ERROR+0.75*ERROR)*1)+Math.log((0.75*NOT_ERROR+0.25*ERROR)*1),
//				0+Math.log((0.25*NOT_ERROR+0.75*ERROR)*1)+Math.log((0.75*NOT_ERROR+0.25*ERROR)*1)
//			};

	}



	@Test
	public void testCalculateLikelihoodCustomSpectrum() {
		String[] seqs = new String[]{
				"AAAC",
				"AACT",
				"ACGT"
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 1);
		Spectrum spectrum = spectrumModel.getSpectrum(0);
		for (int i = 0; i < spectrum.getLength(); i++) {
			double[] freqs = new double[]{1-(0.1*i*3), 0.1*i, 0.1*i, 0.1*i};
			spectrum.resetFrequencies(i, freqs);
			System.out.println("SITE: "+i +"\t"+  Arrays.toString(spectrum.getFrequencies(i)));
		}
//		spectrumModel.setSpectrum(0, spectrum);
		
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double[] eachLikelihood = likelihood.getEachLikelihood();
		System.out.println(Arrays.toString(eachLikelihood));
		//Site1: 1, 0, 0, 0
		//Site2: 0.7, 0.1, 0.1, 0.1
		//Site3: 0.4, 0.2, 0.2, 0.2
		//Site4: 0.1, 0.3, 0.3, 0.3
		double[] expecteds = new double[]{ 
				//MMMD
				Math.log(1*NOT_ERROR+0*ERROR)+Math.log(0.7*NOT_ERROR+0.3*ERROR)+
				Math.log(0.4*NOT_ERROR+0.6*ERROR)+Math.log(0.3*NOT_ERROR+0.7*ERROR),
				
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
