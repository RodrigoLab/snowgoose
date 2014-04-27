package test.srp.evolution.spectrum;

import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import srp.evolution.spectrum.AbstractSpectra;
import srp.evolution.spectrum.SpectraParameter;
import srp.evolution.spectrum.SpectraParameter.SpectraType;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter.DefaultBounds;

public class AbstractSpectraTest {

	@Rule
	public ExpectedException exception = ExpectedException.none();

	private AbstractSpectra spectra;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		spectra = new SpectraParameter(SpectraType.RANDOM);
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testExceptionAddBound() throws Exception {
		Bounds<Double> boundary = new DefaultBounds(1.0, 0.0, 4);
		exception.expect(IllegalArgumentException.class);
		spectra.addBounds(boundary);
		
//		
//		Bounds<Double> boundary2 = new DefaultBounds(10.0, 0.0, 4);
//		exception.expect(IllegalArgumentException.class);
//		spectra.addBounds(boundary2);
	}

	
	@Test
	public void testExceptionGetBound() throws Exception {
		
		Bounds<Double> bounds = spectra.getBounds();
		assertEquals(0.0, bounds.getLowerLimit(0), 0);
		assertEquals(1.0, bounds.getUpperLimit(1), 0);
		assertEquals(4, bounds.getBoundsDimension());

		
//		try {

//			exception.expect(NullPointerException.class);
//			spectra.getBounds();
//		} catch (Exception e) {
//			
//		}

	}
	
	@Test
	public void testExceptionGetParameterValue() throws Exception {
		exception.expect(IllegalArgumentException.class);
		spectra.getParameterValue(0);
	}
	
	@Test
	public void testExceptionGetParameterValues() throws Exception {
		exception.expect(IllegalArgumentException.class);
		spectra.getParameterValues();
	}
	
	@Test
	public void testExceptionSetParameterValue() throws Exception {
		exception.expect(IllegalArgumentException.class);
		spectra.setParameterValue(0,0);
	}
	
	@Test
	public void testExceptionSetPValueQuietly() throws Exception {
		exception.expect(IllegalArgumentException.class);
		spectra.setParameterValueQuietly(0, 0);
		
	}
	@Test
	public void testExceptionSetPValueAll() throws Exception {
		exception.expect(IllegalArgumentException.class);
		spectra.setParameterValueNotifyChangedAll(0, 0);
		
	}
	@Test
	public void testExceptionSetDim() throws Exception {
		exception.expect(IllegalArgumentException.class);
		spectra.setDimension(0);
		
	}
	@Test
	public void testExceptionAddDim() throws Exception {
		exception.expect(IllegalArgumentException.class);
		spectra.addDimension(0, 0);
		
	}
	@Test
	public void testExceptionRemoveDim() throws Exception {
		exception.expect(IllegalArgumentException.class);
		spectra.removeDimension(0);
		
	}

	@Test
	public void testExceptionadoptValues() throws Exception {
		exception.expect(IllegalArgumentException.class);
		spectra.adoptValues(null);
		
	}
	
}

	