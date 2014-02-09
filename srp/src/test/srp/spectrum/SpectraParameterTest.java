package test.srp.spectrum;


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

import srp.spectrum.SpectraParameter;

public class SpectraParameterTest {
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
	public void testConstructor()  {
		double[] frequencies = new double[]{0.25, 0.25, 0.25};
		Assert.assertNotNull(frequencies);
		exception.expect(IllegalArgumentException.class);
		exception.expectMessage("Frequencies should have 4 elements, frequencies.length= 3");
		SpectraParameter spectra = new SpectraParameter(frequencies);
//		Assert.fail();
	}
	
	
//	@Test(expected=IllegalArgumentException.class)
	@Test
	public void testConstructor1()  {		
		double[] frequencies = new double[]{0.2, 0.2, 0.2, 0.2, 0.2};
		exception.expect(IllegalArgumentException.class);
		exception.expectMessage("Frequencies should have 4 elements, frequencies.length= 5");
		new SpectraParameter(frequencies);
	}
	
	@Test
	public void testConstructor2()  {
		double[] frequencies = new double[]{0.25, 0.25, 0.25, 0.15};
		exception.expect(IllegalArgumentException.class);
		exception.expectMessage("Frequencies do not sum to 1, they sum to 0.9");
		new SpectraParameter(frequencies);
	}
	@Test
	public void testConstructor3()  {
		double[] frequencies = new double[]{0.25, 0.25, 0.25, 0.35};
		exception.expect(IllegalArgumentException.class);
		exception.expectMessage("Frequencies do not sum to 1, they sum to 1.1");
		new SpectraParameter(frequencies);
	}
	@Test
	public void testConstructor4()  {
		double[] frequencies = new double[]{0.25, -0.25, 0.5, 0.5};
		exception.expect(IllegalArgumentException.class);
		exception.expectMessage("Frequencies out of bounds 0 < f < 1\t"+ Arrays.toString(frequencies));
		new SpectraParameter(frequencies);

	}
	@Test
	public void testConstructor5()  {
		double[] frequencies = new double[]{1.1, -0.3, 0.1, 0.1};
		exception.expect(IllegalArgumentException.class);
		exception.expectMessage("Frequencies out of bounds 0 < f < 1\t"+ Arrays.toString(frequencies));
		new SpectraParameter(frequencies);
	}	
	
	@Test
	public void testGetSetFrequencies() throws Exception {
		SpectraParameter spectra = new SpectraParameter(0);
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
	
		
	
}
