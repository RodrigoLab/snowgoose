package test.likelihood;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import likelihood.ZOLD_ShortReadLikelihood;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import core.DataImporter;
import dr.evolution.alignment.Alignment;
import dr.evolution.sequence.Sequences;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;

public class ZOLD_ShortReadLikelihoodTest {

	private static final int BINOMIAL = ZOLD_ShortReadLikelihood.BINOMIAL;
	private static final int POISSON = ZOLD_ShortReadLikelihood.POISSON;
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
		haplotypes.add("AACCGGTT"); // 0 mismatch
		ZOLD_ShortReadLikelihood srL = new ZOLD_ShortReadLikelihood(shortReads, haplotypes);
		double logLikelihood = srL.getLogLikelihoodSelect(POISSON);
		double expected =  -0.0856; //dpois(0,8*E,log=T)
		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);
		
		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		expected = -0.086061253223681313806; //dbinom(0,8,E,log=T)
		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);
		
		haplotypes.clear();
		haplotypes.add("TTGGCCAA"); // 8 mismatches
		srL.setHaplotypes(haplotypes);
		logLikelihood = srL.getLogLikelihoodSelect(POISSON);
		expected = -30.3547628694208; //dpois(8,8*E,log=T)
		assertEquals("8 mismatches",expected, logLikelihood, 1e-10);
		
		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		expected = -36.300092300114215504; //dbinom(0,8,E,log=T)
		assertEquals("8 mismatches",expected, logLikelihood, 1e-10);
	}
	
	@Test
	public void testCalculateLikelihood2(){

		shortReads.add("AAAACCCC");
		haplotypes.add("AAAACCCCAAAA"); //0 2 4 6 8 mismatches
		ZOLD_ShortReadLikelihood srL = new ZOLD_ShortReadLikelihood(shortReads, haplotypes);
		
		double logLikelihood = srL.getLogLikelihoodSelect(POISSON);
		double expected = -0.0819407854607265; // log(sum(dpois(seq(0,8,by=2),8*E)))
		
		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		expected = -0.082790211573434788206; // log(sum(dbinom(seq(0,8,by=2),8,E)))
		assertEquals(expected, logLikelihood, 1e-10);
		
		
		haplotypes.clear();
		haplotypes.add("AAAACCCCCCCC"); // 0 1 2 3 4 mismatches
		srL.setHaplotypes(haplotypes);
		logLikelihood = srL.getLogLikelihoodSelect(POISSON);
		expected = -3.56647868713894e-08; // log(sum(dpois(0:4,8*E)))
		assertEquals(expected, logLikelihood, 1e-10);

		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		expected = -7.6461084423155947821e-09; //  log(sum(dbinom(0:4,8,E)))
		assertEquals(expected, logLikelihood, 1e-10);

		
	}

	@Test
	public void testCalculateLikelihood3() {
		
		
		shortReads.add("AAAACCCC");
		haplotypes.add("AAAACCCCAAAA"); //0 2 4 6 8 mismatches
		haplotypes.add("AAAACCCCCCCC"); // 0 1 2 3 4 mismatches
		ZOLD_ShortReadLikelihood srL = new ZOLD_ShortReadLikelihood(shortReads, haplotypes);
		
		double logLikelihood = srL.getLogLikelihoodSelect(POISSON);
		double expected = 0.653015821111653; // log(sum(dpois(c(seq(0,8,by=2),0:4),8*E)))
		assertEquals("multiple haplotypes", expected, logLikelihood, 1e-10);
		
		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		expected =  0.65260860360614092457; // log(sum(dbinom(c(seq(0,8,by=2),0:4),8,E)))
		assertEquals("multiple haplotypes", expected, logLikelihood, 1e-10);
		
		
		shortReads.add("AACCCAA");	//		4 3 1 1 3 4 mismatches
		 							// and	4 3 2 3 4 4 mismatches
		srL.setShortRead(shortReads);
		logLikelihood = srL.getLogLikelihoodSelect(POISSON);
		expected  =  -1.299909707366587952; 
		// log(sum(dpois(c(seq(0,8,by=2),0:4),8*E))) + log(sum(dpois(c(4,3,1,1,3,4,4,3,2,3,4,4),7*E)))
		assertEquals("multiple short reads", expected, logLikelihood, 1e-10);
		
		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		expected =  -1.2931321407960054692;
		// log(sum(dbinom(c(seq(0,8,by=2),0:4),8,E))) + log(sum(dbinom(c(4,3,1,1,3,4,4,3,2,3,4,4),7,E)))
		assertEquals("multiple haplotypes", expected, logLikelihood, 1e-10);
		
		haplotypes.add("TTTTGGGGATTGGACAC");
//		                          8,8,8,8,8,7,7,6,6,6	mismatches
//		                          7,7,6,6,7,7,7,5,5,6,5 mismatches
		srL.setHaplotypes(haplotypes);
		logLikelihood = srL.getLogLikelihoodSelect(POISSON);
		expected  = -1.2999093162319743655; //log(sum(dpois(c(seq(0,8,by=2),0:4, 8,8,8,8,8,7,7,6,6,6),8*E)))+ log(sum(dpois(c(4,3,1,1,3,4,4,3,2,3,4,4,  7,7,6,6,7,7,7,5,5,6,5),7*E)))
		assertEquals(expected, logLikelihood, 1e-10);

		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		expected =  -1.2931320799871470761;
		// log(sum(dpois(c(seq(0,8,by=2),0:4, 8,8,8,8,8,7,7,6,6,6),8*E)))+ log(sum(dpois(c(4,3,1,1,3,4,4,3,2,3,4,4,  7,7,6,6,7,7,7,5,5,6,5),7*E)))
		assertEquals("multiple haplotypes", expected, logLikelihood, 1e-10);

		
	}
	
	@Test
	public void testCompoundLikelihood() {

		

		shortReads.add("AAAACCCC");
		shortReads.add("AACCCAA");
		
		haplotypes.add("AAAACCCCAAAA"); //0 2 4 6 8 mismatches
		haplotypes.add("AAAACCCCCCCC"); // 0 1 2 3 4 mismatches
		haplotypes.add("TTTTGGGGATTGGACAC");
		
		ZOLD_ShortReadLikelihood srL = new ZOLD_ShortReadLikelihood(shortReads, haplotypes);
		
		double logLikelihood = srL.getLogLikelihood();
		double expected =  -1.2999093162319743655; //log(sum(dpois(c(seq(0,8,by=2),0:4, 8,8,8,8,8,7,7,6,6,6),8*E)))+ log(sum(dpois(c(4,3,1,1,3,4,4,3,2,3,4,4,  7,7,6,6,7,7,7,5,5,6,5),7*E)))
		assertEquals("multiple short reads", expected, logLikelihood, 1e-10);
		
        List<Likelihood> likelihoods = new ArrayList<Likelihood>();        
        likelihoods.add(srL);
        Likelihood CSRLikelihood = new CompoundLikelihood(0, likelihoods);
        CSRLikelihood.setId(ZOLD_ShortReadLikelihood.NAME);
//        System.out.println(((CompoundLikelihood) CSRLikelihood).getDiagnosis());
//        System.out.println(((CompoundLikelihood) CSRLikelihood).getReport());
//        System.out.println(CSRLikelihood.getId());
//        System.out.println(CSRLikelihood.prettyName());
//        System.out.println(CSRLikelihood.getModel().getModelName());
        logLikelihood = CSRLikelihood.getLogLikelihood();
        assertEquals("multiple short reads", expected, logLikelihood, 1e-10);
        

	}
	
	@Test
	public void testDifferentImplementation100Reads(){
		
		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/Stage0/";
		
		if( Files.exists(Paths.get(dataDir)) ){
	
			String trueAlignmentFile = "121101_true_seqs.fasta";
			String shortReadFile = "121101_short_reads_100.fasta";
			
			DataImporter dataImporter = new DataImporter(dataDir);
			Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
			Sequences shortReads = dataImporter.importSequence(shortReadFile);
			
			ZOLD_ShortReadLikelihood srl = new ZOLD_ShortReadLikelihood(shortReads, trueAlignment);
		
			double expected = -106.642313524702; //100 reads
			
			// calculateShoreReadLikelihood1 ~ 14500-15000
			// calculateShoreReadLikelihood2 ~ 800-900
			
			double actual = srl.getLogLikelihoodSelect(0);
			assertEquals("Compare to python version, 100 reads", expected, actual, 1e-10);
			actual = srl.getLogLikelihoodSelect(1);
			assertEquals("Compare to python version, 100 reads", expected, actual, 1e-10);
//			actual = srl.getLogLikelihoodSelect(2);
//			assertEquals("Compare to python version, 100 reads", expected, actual, 1e-10);
		
			
		} else {
			fail("change hardcoded dataDir: " + dataDir);
		}
	}

	

	@Test
	public void testDifferentImplementation10Reads(){
		
		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/Stage0/";
		
		if( Files.exists(Paths.get(dataDir)) ){
	
			String trueAlignmentFile = "121101_true_seqs.fasta";
			String shortReadFile = "121101_short_reads_10.fasta";
			
			DataImporter dataImporter = new DataImporter(dataDir);
			Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
			Sequences shortReads = dataImporter.importSequence(shortReadFile);
			
			ZOLD_ShortReadLikelihood srl = new ZOLD_ShortReadLikelihood(shortReads, trueAlignment);
		
			double expected = -11.85523370076271; //10  reads
			
			double actual = srl.getLogLikelihoodSelect(0);
			assertEquals("Compare to python version, 10 reads", expected, actual, 1e-10);
			actual = srl.getLogLikelihoodSelect(1);
			assertEquals("Compare to python version, 10 reads", expected, actual, 1e-10);
//			actual = srl.getLogLikelihoodSelect(2);
//			assertEquals("Compare to python version, 100 reads", expected, actual, 1e-10);
		
			
		} else {
			fail("change hardcoded dataDir: " + dataDir);
		}
	}

}
