package test.srp.spectrum;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import srp.evolution.spectrum.CategorySpectraParameter;
import srp.evolution.spectrum.CategorySpectraParameter.CategoryType;

public class CategorySpectraParameterTest {

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
	public void testConstructorException() throws Exception {
		exception.expect(IllegalArgumentException.class);
		CategorySpectraParameter spectra = new CategorySpectraParameter(100);
		
	}

	@Test
	public void testConstructorException2() throws Exception {
		exception.expect(IllegalArgumentException.class);
		CategorySpectraParameter spectra = new CategorySpectraParameter(-1);
		
	}
	
	@Test
	public void testConstructor() throws Exception {
		
		CategorySpectraParameter spectra = new CategorySpectraParameter(0);
		int cat = spectra.getCategory();
		int expected = 0;
		assertEquals(expected, cat);
		
		double freq;
		double[] expecteds = new double[]{0.97, 0.01, 0.01, 0.01};
		double[] frequencies = spectra.getFrequencies();
		assertArrayEquals(expecteds, frequencies, 0);
		for (int i = 0; i < expecteds.length; i++) {
			freq = spectra.getFrequency(i);
			assertEquals(expecteds[i], freq, 0);
		}
		
		
		spectra = new CategorySpectraParameter(4);
		cat = spectra.getCategory();
		expected = 4;
		assertEquals(expected, cat);
		
		expecteds = new double[]{0.49, 0.49, 0.01, 0.01};
		frequencies = spectra.getFrequencies();
		assertArrayEquals(expecteds, frequencies, 0);
		for (int i = 0; i < expecteds.length; i++) {
			freq = spectra.getFrequency(i);
			assertEquals(expecteds[i], freq, 0);
		}
	}
	
	
	@Test
	public void testConstructor2() throws Exception {
		CategorySpectraParameter spectra;
		int[] count = new int[4];
		for (int i = 0; i < 100; i++) {
			spectra = new CategorySpectraParameter(CategoryType.SINGLE);
			int cat = spectra.getCategory();
			count[cat]++;
			assertTrue("Category outside range: "+cat, cat >= 0 && cat <= 3);
		}
		for (int i = 0; i < count.length; i++) {
			assertTrue("i="+i+", Count[i]="+count[i], count[i]>0);
		}
		
		count = new int[10];
		for (int i = 0; i < 100; i++) {
			spectra = new CategorySpectraParameter(CategoryType.TWOWAYS);
			int cat = spectra.getCategory();
			count[cat]++;
			assertTrue("Category outside range: "+cat, cat >= 0 && cat <= 9);
		}
		for (int i = 0; i < count.length; i++) {
			assertTrue("i="+i+", Count[i]="+count[i], count[i]>0);
		}
	}
	
	
	
	@Test
	public void testSetCategory() throws Exception {
		CategorySpectraParameter spectra = new CategorySpectraParameter(0);
		int cat = spectra.getCategory();
		int expected = 0;
		assertEquals(expected, cat);
		
		spectra.setCategory(1);
		
		double freq;
		double[] expecteds = new double[]{0.01, 0.97, 0.01, 0.01};
		double[] frequencies = spectra.getFrequencies();
		assertArrayEquals(expecteds, frequencies, 0);
		for (int i = 0; i < expecteds.length; i++) {
			freq = spectra.getFrequency(i);
			assertEquals(expecteds[i], freq, 0);
		}
		
		
		spectra = new CategorySpectraParameter(5);
		cat = spectra.getCategory();
		expected = 5;
		assertEquals(expected, cat);
		
		expecteds = new double[]{0.49, 0.01, 0.49, 0.01};
		frequencies = spectra.getFrequencies();
		assertArrayEquals(expecteds, frequencies, 0);
		for (int i = 0; i < expecteds.length; i++) {
			freq = spectra.getFrequency(i);
			assertEquals(expecteds[i], freq, 0);
		}
	}
	
	
	@Test
	public void testStoreRestore() throws Exception {
		//TODO implement testStoreRestore()
	}
	
	
	
	/**
CategorySpectraParameter()
CategorySpectraParameter(CategoryType)
CategorySpectraParameter(int)
setCategory(int)
getCategory()
getFrequencies()
getFrequency(int)


	 */
}
