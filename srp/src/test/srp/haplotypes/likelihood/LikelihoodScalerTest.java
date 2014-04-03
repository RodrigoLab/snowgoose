package test.srp.haplotypes.likelihood;

import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.likelihood.haplotypes.LikelihoodScaler;

public class LikelihoodScalerTest {

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
	public void testGetLogLikelihood() {
		final double C = 1e-200;
		final double LOG_C = Math.log(C);
		
		LikelihoodScaler liS = new LikelihoodScaler(LOG_C);
		
		double expected = 0;
		for (int i = 0; i < 100; i++) {
			double data = Math.random();
			expected += data;
			liS.addLogProb(Math.log(data));
		}
		expected = Math.log(expected);
		
		double actual = liS.getLogLikelihood();
		assertEquals(expected, actual, 1e-10);
		
		
	}

	@Test
	public void testGetLogLikelihoodAdd() {
		final double C = 1e-200;
		final double LOG_C = Math.log(C);
		
		LikelihoodScaler liS = new LikelihoodScaler(LOG_C);
		
		double expected = 0;
		for (int i = 0; i < 100; i++) {
			double data = Math.random();
			expected += data;
			double logData = liS.scale(Math.log(data));
			liS.add(logData);
		}
		expected = Math.log(expected);
		
		double actual = liS.getLogLikelihood();
		assertEquals(expected, actual, 1e-10);
		
		
	}
	
	
}
