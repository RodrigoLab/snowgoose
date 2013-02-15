package test.likelihood;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import likelihood.ShortReadLikelihood;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import alignment.AlignmentMapping;
import alignment.AlignmentMatrix;
import alignment.AlignmentUtils;

import core.DataImporter;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.Sequences;
import dr.evolution.util.Taxon;
import dr.inference.mcmc.MCMCCriterion;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;

public class ShortReadLikelihoodTest {

	private static final int BINOMIAL = ShortReadLikelihood.BINOMIAL;

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
	public void testCalculateLikelihoodIdentical() {

		String[] seqs = new String[]{"AACCGGTT"};
	
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		SimpleAlignment hap = AlignmentUtils.createAlignment(seqs);
		
		AlignmentMatrix aMatrix = new AlignmentMatrix(aMap, hap);
		ShortReadLikelihood srL = new ShortReadLikelihood(aMap, aMatrix);
		
		double logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		double expected = -0.086061253223681313806; //dbinom(0,8,E,log=T)
		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);
		
		seqs = new String[]{
				".....AACCGGTT.",
				".....AACCGGTT...",
				".....AACCGGTT.....",
				".....AACCGGTT.......",
				};
		aMap = AlignmentUtils.createAlignmentMapping(seqs);
		seqs = new String[]{"XXXXXAACCGGTTYYYYYYYY"};
		hap = AlignmentUtils.createAlignment(seqs);
		
		aMatrix = new AlignmentMatrix(aMap, hap);
		srL = new ShortReadLikelihood(aMap, aMatrix);
		
		
		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		expected = -0.086061253223681313806*4; //dbinom(0,8,E,log=T)
		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);
		
	}

	@Test
	public void testCalculateLikelihoodSomething() {
		String[] seqs = new String[]{
				"AAAAAAAAATGTGTTTT....",
				".....CCCCCCCCCCCCCCCCCCCTTTTCCCC....",
				"..........GGGGGGGGGGGGGGCGCGTATAGGGG",
				"...............TTTTTTTTTACACTATA....",
//				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG"
//				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		String[] haps = new String[]{
//				"AAAAACCCCCGGGGGTTTTTACGTACACTATATATA"
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		AlignmentMatrix aMatrix = new AlignmentMatrix(aMap, hapAlignment);

		ShortReadLikelihood srL = new ShortReadLikelihood(aMap, aMatrix);
		double logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		
		System.out.println(aMatrix.getHaplotype(0) );
		System.out.println("init:\t"+logLikelihood);
		System.out.println();
		
//		aMatrix = new AlignmentMatrix(aMap, 1);
//		srL = new ShortReadLikelihood(aMap, aMatrix);
//		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
//		System.out.println("Random init:\t"+logLikelihood);
		
		
//		aMatrix.swapSrp(0, 0, 10, 1);
//		System.out.println(aMatrix.getHaplotype(0) );
//		srL = new ShortReadLikelihood(aMap, aMatrix);
//		System.out.println("Likelihood:\t"+srL.getLogLikelihood()+"\n");
		aMatrix = new AlignmentMatrix(aMap, 5);
		for (int i = 0; i < aMatrix.getNoHap(); i++) {
			aMatrix.randomSeq(i);
			System.out.println(aMatrix.getHaplotype(i) );
		}

		srL = new ShortReadLikelihood(aMap, aMatrix);
		System.out.println("Likelihood:\t"+srL.getLogLikelihood()+"\n");
		
		double likelihood = srL.getLogLikelihood(); 
		
		MCMCCriterion criterion = new MCMCCriterion();
        double hastingsRatio = 0.0;
        double[] logr = {-Double.MAX_VALUE};
        
		for (int i = 0; i < 1e5; i++) {
			
			aMatrix.swapBase();
			srL.updatehaplotypesChars(aMatrix);
//			System.out.print(srL.getLogLikelihood() +"\t" );
//			srL = new ShortReadLikelihood(aMap, aMatrix);
//			System.out.println(srL.getLogLikelihood());
			double newL = srL.getLogLikelihood();
			
			
			
			boolean accept = criterion.accept(likelihood, newL, hastingsRatio, logr);
//			System.out.println(likelihood +"\t"+newL +"\t"+ logr[0] +"\t"+  accept);
			if(accept){
				likelihood = newL;
			}
			else{
				aMatrix.reject();
				
			}
			
				
		}
		for (int i = 0; i < aMatrix.getNoHap(); i++) {
			System.out.println(aMatrix.getHaplotype(i) );
		}
		srL = new ShortReadLikelihood(aMap, aMatrix);
		System.out.println("Likelihood:\t"+srL.getLogLikelihood()+"\n");
		

//AAAAACCCCCGGGCCTTCGCGTCCACTTTATAGGGG
//Likelihood:	-142.41751922532444


//		randomSeq

	}	
	@Test
	public void testCalculateLikelihoodAllDiff() {

		String[] seqs = new String[]{"AACCGGTT"};
		SimpleAlignment srp = AlignmentUtils.createAlignment(seqs);
		AlignmentMapping aMap = new AlignmentMapping(srp);
		
		seqs = new String[]{"TTGGCCAA"};
		SimpleAlignment hap = AlignmentUtils.createAlignment(seqs);
		
		AlignmentMatrix aMatrix = new AlignmentMatrix(aMap, hap);
		ShortReadLikelihood srL = new ShortReadLikelihood(aMap, aMatrix);
		
		
		double logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		double expected = -36.300092300114215504; //dbinom(0,8,E,log=T)
		assertEquals("8 mismatches",expected, logLikelihood, 1e-10);
		
//		
//		srp.removeSequence(0);
//		srp.addSequence(new Sequence(t, ".....TTGGCCAA......."));
//		aMap = new AlignmentMapping(srp);
//		System.out.println(aMap.getFullSrp(0));
//		System.out.println(aMap.getFragmentSrp(0));
//		srL = new ShortReadLikelihood(aMap, hap);
//		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
//		assertEquals("8 mismatches",expected, logLikelihood, 1e-10);

	}
	
	@Test
	public void testCalculateLikelihood2(){

		SimpleAlignment srp = new SimpleAlignment();
		Taxon t = new Taxon("TEST");
		Sequence s = new Sequence(t, "AACCGGTT");
		srp.addSequence(s);		
		AlignmentMapping aMap = new AlignmentMapping(srp);
		
		SimpleAlignment hap = new SimpleAlignment();
		hap.addSequence(new Sequence(t, "TTGGCCAA"));
		
		ShortReadLikelihood srL = new ShortReadLikelihood(aMap, srp);
		
		double logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		double expected = -0.086061253223681313806; //dbinom(0,8,E,log=T)
		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);
		

		
		
		shortReads.add("AAAACCCC");
		haplotypes.add("AAAACCCCAAAA"); //0 2 4 6 8 mismatches
		ShortReadLikelihood srL = new ShortReadLikelihood(shortReads, haplotypes);
		
		double logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		double expected = -0.0819407854607265; // log(sum(dpois(seq(0,8,by=2),8*E)))
		
		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		expected = -0.082790211573434788206; // log(sum(dbinom(seq(0,8,by=2),8,E)))
		assertEquals(expected, logLikelihood, 1e-10);
		
		
		haplotypes.clear();
		haplotypes.add("AAAACCCCCCCC"); // 0 1 2 3 4 mismatches
		srL.setHaplotypes(haplotypes);
		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
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
		ShortReadLikelihood srL = new ShortReadLikelihood(shortReads, haplotypes);
		
		double logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		double expected = 0.653015821111653; // log(sum(dpois(c(seq(0,8,by=2),0:4),8*E)))
		assertEquals("multiple haplotypes", expected, logLikelihood, 1e-10);
		
		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
		expected =  0.65260860360614092457; // log(sum(dbinom(c(seq(0,8,by=2),0:4),8,E)))
		assertEquals("multiple haplotypes", expected, logLikelihood, 1e-10);
		
		
		shortReads.add("AACCCAA");	//		4 3 1 1 3 4 mismatches
		 							// and	4 3 2 3 4 4 mismatches
		srL.setShortRead(shortReads);
		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
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
		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
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
		
		ShortReadLikelihood srL = new ShortReadLikelihood(shortReads, haplotypes);
		
		double logLikelihood = srL.getLogLikelihood();
		double expected =  -1.2999093162319743655; //log(sum(dpois(c(seq(0,8,by=2),0:4, 8,8,8,8,8,7,7,6,6,6),8*E)))+ log(sum(dpois(c(4,3,1,1,3,4,4,3,2,3,4,4,  7,7,6,6,7,7,7,5,5,6,5),7*E)))
		assertEquals("multiple short reads", expected, logLikelihood, 1e-10);
		
        List<Likelihood> likelihoods = new ArrayList<Likelihood>();        
        likelihoods.add(srL);
        Likelihood CSRLikelihood = new CompoundLikelihood(0, likelihoods);
        CSRLikelihood.setId(ShortReadLikelihood.NAME);
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
			
			ShortReadLikelihood srl = new ShortReadLikelihood(shortReads, trueAlignment);
		
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
			
			ShortReadLikelihood srl = new ShortReadLikelihood(shortReads, trueAlignment);
		
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
