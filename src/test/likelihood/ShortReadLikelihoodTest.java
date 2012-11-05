package test.likelihood;

import static org.junit.Assert.*;

import java.util.ArrayList;

import likelihood.ShortReadLikelihood;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

public class ShortReadLikelihoodTest {

	ArrayList<String> shortReads;
	ArrayList<String> haplotypes;
	

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {

	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {

		shortReads = new ArrayList<>();
		haplotypes = new ArrayList<>();
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testCalculateLikelihood() {
		
		shortReads.add("AACCGGTT");
		haplotypes.add("AACCGGTT"); // 0 mis-matches
		ShortReadLikelihood srL = new ShortReadLikelihood(shortReads, haplotypes);
		double logLikelihood = srL.calculateLikelihood();
		double expected =  -0.0856; //dpois(0,8*E,log=T)
		assertEquals(expected, logLikelihood, 1e-10);
	
		
		haplotypes.clear();
		haplotypes.add("TTGGCCAA"); // 8 mis-match
		srL.setHaplotypes(haplotypes);
		logLikelihood = srL.calculateLikelihood();
		expected = -30.3547628694208; //dpois(8,8*E,log=T)
		assertEquals(expected, logLikelihood, 1e-10);
	}
	@Test
	public void testCalculateLikelihood2(){

		shortReads.add("AAAACCCC");
		haplotypes.add("AAAACCCCAAAA"); //0 2 4 6 8 mismatches
		ShortReadLikelihood srL = new ShortReadLikelihood(shortReads, haplotypes);
		
		double logLikelihood = srL.calculateLikelihood();
		double expected = -0.0819407854607265; // log(sum(dpois(seq(0,8,by=2),8*E)))
		assertEquals(expected, logLikelihood, 1e-10);

		haplotypes.clear();
		haplotypes.add("AAAACCCCCCCC"); // 0 1 2 3 4 mismatches
		srL.setHaplotypes(haplotypes);
		logLikelihood = srL.calculateLikelihood();
		expected = -3.56647868713894e-08; // log(sum(dpois(0:4,8*E)))
		assertEquals(expected, logLikelihood, 1e-10);

		
	}

	@Test
	public void testCalculateLikelihood3() {
		
		
		shortReads.add("AAAACCCC");
		haplotypes.add("AAAACCCCAAAA"); //0 2 4 6 8 mismatches
		haplotypes.add("AAAACCCCCCCC"); // 0 1 2 3 4 mismatches
		ShortReadLikelihood srL = new ShortReadLikelihood(shortReads, haplotypes);
		
		double logLikelihood = srL.calculateLikelihood();
		double expected = 0.653015821111653; // log(sum(dpois(c(seq(0,8,by=2),0:4),8*E)))
		assertEquals(expected, logLikelihood, 1e-10);
		
		shortReads.add("AACCCAA");	//		4 3 1 1 3 4 mismatches
		 							// and	4 3 2 3 4 4 mismatches
		logLikelihood = srL.calculateLikelihood();
		expected += -1.9529255284782407465; // log(sum(dpois(c(4,3,1,1,3,4,4,3,2,3,4,4),7*E)))
		assertEquals(expected, logLikelihood, 1e-10);
		
		
		haplotypes.add("TTTTGGGGATTGGACAC");
//		                          8,8,8,8,8,7,7,6,6,6	mismatches
//		                          7,7,6,6,7,7,7,5,5,6,5 mismatches
		logLikelihood = srL.calculateLikelihood();
		expected  = -1.2999093162319743655; //log(sum(dpois(c(seq(0,8,by=2),0:4, 8,8,8,8,8,7,7,6,6,6),8*E)))+ log(sum(dpois(c(4,3,1,1,3,4,4,3,2,3,4,4,  7,7,6,6,7,7,7,5,5,6,5),7*E)))
		assertEquals(expected, logLikelihood, 1e-10);

		
	}
}
