package test.likelihood;

import static org.junit.Assert.*;
import likelihood.LikelihoodScaler;
import likelihood.LikelihoodUtils;

import org.apache.commons.math3.stat.StatUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

public class LikelihoodUtilsTest {

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
	public void testHamDist() {
		String s1 = "CATGTGCAGCAGCAGGCGCGACTCTGCAGTCGTACTGCTGACTCTGACTCATCTACTACTACG";
		String s2 = "CATGTGCAGCAGCAGGCGCGACTCTGCAGTCGTACTGCTGACTCTGACTCATCTACTACCTAC";
		int actual = LikelihoodUtils.hamDist(s1, s2);
		assertEquals(4, actual);
		
		s2 = "CATGTGCAGCAGCAGGCGCGACTCTGCAGTCGTACTGCTGACTCTGACTCATCTAGGGGCTAC";
		actual = LikelihoodUtils.hamDist(s1, s2);
		assertEquals(8, actual);
		
	}

}
