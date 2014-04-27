package test.srp.evolution.spectrum;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.evolution.spectrum.CategorySpectrum;
import srp.evolution.spectrum.CategorySpectrumAlignmentModel;

public class CategorySpectrumAlignmentModelTest {

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
		CategorySpectrumAlignmentModel alignment = new CategorySpectrumAlignmentModel(10, 3);
		assertEquals(10, alignment.getSpectrumLength());
		assertEquals(3, alignment.getSpectrumCount());
		
		
		double[][] expectedsFreq= new double[][]{
			{0.97, 0.01, 0.01, 0.01},
			{0.01, 0.97, 0.01, 0.01},
			{0.01, 0.01, 0.97, 0.01},
			{0.01, 0.01, 0.01, 0.97},
		};
		
		for (int i = 0; i < alignment.getSpectrumCount(); i++) {
			CategorySpectrum spectrum = alignment.getSpectrum(i);
			
			for (int j = 0; j < alignment.getSpectrumLength(); j++) {
				int cat = spectrum.getSpectra(j).getCategory();
				double[] freq = alignment.getSpecturmFrequencies(i, j);
					assertArrayEquals(expectedsFreq[cat], freq, 0);
			}
		}
		
		
	}
	/**
	 * 
	 addSpectrum(CategorySpectrum)
getSpectrumCount()
getSpectrum(int)
getSpectrumString(int)
removeSpectrum(int)
	 */
	
}
