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

	ArrayList<String> shortReads = new ArrayList<>();
	ArrayList<String> haplotypes = new ArrayList<>();

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {

	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {

		shortReads.add("AACCGGTT");
		shortReads.add("ACTGG");

		haplotypes.add("CGATGTGTTCTTGGAATCACCTACCCTAGTGAC");
		haplotypes.add("AAAAAACACACACAACCGGTTACGTGTGTGTGT");
		haplotypes.add("AAAAAACACACACTGCATTGGCCAAGTGTGTGT");
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testCalculateLikelihood() {
		ShortReadLikelihood srL = new ShortReadLikelihood(shortReads, haplotypes);
		double logLikelihood = srL.calculateLikelihood();
		System.out.println(logLikelihood);
	}

}
