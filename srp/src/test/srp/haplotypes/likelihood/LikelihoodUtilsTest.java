package test.srp.haplotypes.likelihood;

import static org.junit.Assert.*;

import org.apache.commons.math3.stat.StatUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.haplotypes.likelihood.LikelihoodScaler;
import srp.haplotypes.likelihood.LikelihoodUtils;

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
	@Test
	public void testHamDistChar() {
		String s1 = "CATGTGCAGCAGCAGGCGCGACTCTGCAGTCGTACTGCTGACTCTGACTCATCTACTACTACG";
		char[] c1 = s1.toCharArray();
		String s2 = "CATGTGCAGCAGCAGGCGCGACTCTGCAGTCGTACTGCTGACTCTGACTCATCTACTACCTAC";
		char[] c2 = s2.toCharArray();
		int actual = LikelihoodUtils.Dist(0, c1.length, c1, c2);
		assertEquals(4, actual);
		
		c2 = "CATGTGCAGCAGCAGGCGCGACTCTGCAGTCGTACTGCTGACTCTGACTCATCTAGGGGCTAC".toCharArray();
		actual = LikelihoodUtils.Dist(0, c1.length, c1, c2);
		assertEquals(8, actual);
		
	}
}
