package test.srp.evolution.spectrum;


import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import srp.evolution.spectrum.SpectraParameter;
import srp.evolution.spectrum.SpectraParameter.SpectraType;

public class SpectraParameterTest {
	
	
	private static final int DIMENSION = SpectraParameter.DIMENSION;
	@Rule
	public ExpectedException exception = ExpectedException.none();

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
	public void testConstructorNot4Elements()  {
		
		assertEquals(4, DIMENSION);
		double[] frequencies = new double[]{0.25, 0.25, 0.25};
		Assert.assertNotNull(frequencies);
		exception.expect(IllegalArgumentException.class);
		exception.expectMessage("Frequencies should have 4 elements, frequencies.length= 3");
		SpectraParameter spectra = new SpectraParameter(frequencies);
//		Assert.fail();
	}
	
	
//	@Test(expected=IllegalArgumentException.class)
	@Test
	public void testConstructorNot4Elements2()  {		
		double[] frequencies = new double[]{0.2, 0.2, 0.2, 0.2, 0.2};
		exception.expect(IllegalArgumentException.class);
		exception.expectMessage("Frequencies should have 4 elements, frequencies.length= 5");
		new SpectraParameter(frequencies);
	}
	
	@Test
	public void testConstructorSumLeOne()  {
		double[] frequencies = new double[]{0.25, 0.25, 0.25, 0.15};
		exception.expect(IllegalArgumentException.class);
		exception.expectMessage("Frequencies do not sum to 1, they sum to 0.9");
		new SpectraParameter(frequencies);
	}
	@Test
	public void testConstructorSumGtOne()  {
		double[] frequencies = new double[]{0.25, 0.25, 0.25, 0.35};
		exception.expect(IllegalArgumentException.class);
		exception.expectMessage("Frequencies do not sum to 1, they sum to 1.1");
		new SpectraParameter(frequencies);
	}
	@Test
	public void testConstructorOutOfBound()  {
		double[] frequencies = new double[]{0.25, -0.25, 0.5, 0.5};
		exception.expect(IllegalArgumentException.class);
		exception.expectMessage("Frequencies out of bounds 0 < f < 1\t"+ Arrays.toString(frequencies));
		new SpectraParameter(frequencies);

	}
	@Test
	public void testConstructorOutOfBounds2()  {
		double[] frequencies = new double[]{1.1, -0.3, 0.1, 0.1};
		exception.expect(IllegalArgumentException.class);
		exception.expectMessage("Frequencies out of bounds 0 < f < 1\t"+ Arrays.toString(frequencies));
		new SpectraParameter(frequencies);
	}	
	
	
	@Test
	public void testSpectraType() throws Exception {
		SpectraParameter spectra = new SpectraParameter(SpectraType.EQUAL);
	}
	
	@Test
	public void testGetSetFrequencies() throws Exception {
		SpectraParameter spectra = new SpectraParameter(SpectraType.EQUAL);
		assertEquals(4, spectra.getDimension());
		for (int i = 0; i < spectra.getDimension(); i++) {
			assertEquals(0.25, spectra.getFrequency(i), 0);
			spectra.setFrequency(i, (i+1)/10.0);
		}
		assertEquals(0.1, spectra.getFrequency(0), 0);
		assertEquals(0.2, spectra.getFrequency(1), 0);
		assertEquals(0.3, spectra.getFrequency(2), 0);
		assertEquals(0.4, spectra.getFrequency(3), 0);

	}
	

	@Test
	public void testConstructorEqual() throws Exception {
		double[] expectedFreq = new double[]{0.25, 0.25, 0.25, 0.25};
		for (int i = 0; i < 100; i++) {
			SpectraParameter spectra = new SpectraParameter(SpectraType.EQUAL);
			assertArrayEquals(expectedFreq , spectra.getFrequencies(), 0);
		}
		
	}
	
	@Test
	public void testConstructorZeroOne() throws Exception {
		for (int i = 0; i < 100; i++) {
			SpectraParameter spectra = new SpectraParameter(
					SpectraType.DOMINANT);
			int count0 = 0;
			int count1 = 0;
			for (int k = 0; k < spectra.getDimension(); k++) {
				double freq = spectra.getFrequency(k);
				if (freq == SpectraParameter.INIT_SMALL) {
					count0++;
				} else if (freq == SpectraParameter.INIT_LARGE) {
					count1++;
				}
			}
			assertEquals(3, count0);
			assertEquals(1, count1);
		}

	}


	@Test
	public void testConstructorRandom() throws Exception {
		for (int i = 0; i < 100; i++) {
			SpectraParameter spectra = new SpectraParameter(SpectraType.RANDOM);
			double[] frequencies = spectra.getFrequencies();
			int count = 0;
			for (int k = 0; k < frequencies.length; k++) {
				for (int k2 = k+1; k2 < frequencies.length; k2++) {
					if (frequencies[k] != frequencies[k2]) {
						count++;
					}
				}

			}
			assertEquals(Arrays.toString(frequencies), 6, count);
		}

	}
	@Test
	public void testGetStoredFrequency() throws Exception {
		
		SpectraParameter spectra = new SpectraParameter(SpectraType.RANDOM);
		double[] freq = new double[]{Math.random(), Math.random(), Math.random(), Math.random()};
		for (int i = 0; i < DIMENSION ; i++) {
			spectra.setFrequency(i, freq[i]);
		}
		for (int i = 0; i < 100; i++) {
			spectra.storeState();
			double[] storedFreq = new double[DIMENSION];
			System.arraycopy(freq, 0, storedFreq, 0, storedFreq.length);
			freq = new double[]{Math.random(), Math.random(), Math.random(), Math.random()};
			for (int f = 0; f < freq.length ; f++) {
				spectra.setFrequency(f, freq[f]);
			}
			assertArrayEquals(freq, spectra.getFrequencies(), 0);
			for (int f = 0; f < storedFreq.length; f++) {
				assertEquals(storedFreq[f], spectra.getStoredFrequency(f), 0);
			}
			
		}		
	}
	
	
	@Test
	public void testSetGetStateLikelihood() throws Exception {
//		fail("Not yet implemented");
	}
	
	@Test
	public void testStoreRestore() throws Exception {
//		fail("Not yet implemented");
	}
	
}
