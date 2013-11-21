package test.srp.spectrum;


import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;

public class SpectrumTest {

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
		Spectrum spectrum = new Spectrum(5);
		assertEquals(5, spectrum.getLength());
		for (int s = 0; s < spectrum.getLength(); s++) {
			for (int p = 0; p < 3; p++) {
				assertEquals(0.25, spectrum.getFrequency(s, p), 0);
			}
		}
	}
	@Test
	public void testConstructor2() throws Exception {
		
		Spectrum spectrum = new Spectrum(5, new double[]{0.1, 0.2, 0.3, 0.4});
		assertEquals(5, spectrum.getLength());
		double[] expectedFreq = new double[]{0.1, 0.2, 0.3, 0.4};
		for (int s = 0; s < spectrum.getLength(); s++) {
			assertArrayEquals(expectedFreq, spectrum.getFrequencies(s), 0);
		}
	}
	
	@Test
		public void testResetFrequencies() throws Exception {
			Spectrum spectrum = new Spectrum(5);
//			SpectraParameter spectra = spectrum.getSpectra(0);
			for (int i = 0; i < 4; i++) {
				spectrum.setFrequencyAt(i, i, (i+1)/10.0);
			}
			assertEquals(0.1, spectrum.getFrequency(0, 0), 0);
			assertEquals(0.2, spectrum.getFrequency(1, 1), 0);
			assertEquals(0.3, spectrum.getFrequency(2, 2), 0);
			assertEquals(0.4, spectrum.getFrequency(3, 3), 0);
			
		}
}
