package test.srp.evolution.spectrum;


import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.evolution.spectrum.Spectrum;

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
//	@Test
//	public void testConstructor2() throws Exception {
//		
//		Spectrum spectrum = new Spectrum(5);
//		assertEquals(5, spectrum.getLength());
//		double[] expectedFreq = new double[]{0.1, 0.2, 0.3, 0.4};
//		for (int s = 0; s < spectrum.getLength(); s++) {
//			spectrum
//			assertArrayEquals(expectedFreq, spectrum.getFrequencies(s), 0);
//		}
//	}

	@Test
	public void testResetSpectra() throws Exception {
		Spectrum spectrum = new Spectrum(5);

		for (int i = 0; i < 4; i++) {
			double[] newFreq = new double[] { (i) / 10.0, (i + 1) / 10.0,
					(i + 2) / 10.0, (i + 3) / 10.0 };
			spectrum.resetSpectra(i, newFreq);
		}
		assertEquals(0.0, spectrum.getFrequency(0, 0), 0);
		assertEquals(0.1, spectrum.getFrequency(0, 1), 0);
		assertEquals(0.2, spectrum.getFrequency(0, 2), 0);
		assertEquals(0.3, spectrum.getFrequency(0, 3), 0);

		assertEquals(0.1, spectrum.getFrequency(1, 0), 0);
		assertEquals(0.2, spectrum.getFrequency(1, 1), 0);
		assertEquals(0.3, spectrum.getFrequency(1, 2), 0);
		assertEquals(0.4, spectrum.getFrequency(1, 3), 0);

		assertEquals(0.2, spectrum.getFrequency(2, 0), 0);
		assertEquals(0.3, spectrum.getFrequency(2, 1), 0);
		assertEquals(0.4, spectrum.getFrequency(2, 2), 0);
		assertEquals(0.5, spectrum.getFrequency(2, 3), 0);

		assertEquals(0.3, spectrum.getFrequency(3, 0), 0);
		assertEquals(0.4, spectrum.getFrequency(3, 1), 0);
		assertEquals(0.5, spectrum.getFrequency(3, 2), 0);
		assertEquals(0.6, spectrum.getFrequency(3, 3), 0);

	}

	@Test
	public void testStoreRestore() throws Exception {
		//TODO implement testStoreRestore()
	}
}
