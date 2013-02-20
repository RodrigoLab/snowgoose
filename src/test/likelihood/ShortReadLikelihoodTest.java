package test.likelihood;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import likelihood.ShortReadLikelihood;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import alignment.AlignmentMapping;
import alignment.Haplotypes;
import alignment.AlignmentUtils;

import core.DataImporter;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.Sequences;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.inference.mcmc.MCMCCriterion;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;

public class ShortReadLikelihoodTest {

	private static final int BINOMIAL = ShortReadLikelihood.BINOMIAL;

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {

	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {

//		shortReads = new ArrayList<>();
//		haplotypes = new ArrayList<>();
	}

	@After
	public void tearDown() throws Exception {
	}

	
	
	@Test
	public void testCalculateLikelihoodIdentical() {

		String[] seqs = new String[]{"AACCGGTT"};
	
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		SimpleAlignment hap = AlignmentUtils.createAlignment(seqs);
		
		Haplotypes aMatrix = new Haplotypes(aMap, hap);
		ShortReadLikelihood srL = new ShortReadLikelihood(aMap, aMatrix);
		
		double logLikelihood = srL.calculateLogLikelihoodSelect(BINOMIAL);
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
		
		aMatrix = new Haplotypes(aMap, hap);
		srL = new ShortReadLikelihood(aMap, aMatrix);
		
		
		logLikelihood = srL.calculateLogLikelihoodSelect(BINOMIAL);
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
//				"AAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		Haplotypes aMatrix = new Haplotypes(aMap, hapAlignment);

		ShortReadLikelihood srL = new ShortReadLikelihood(aMap, aMatrix);
		double logLikelihood = srL.calculateLogLikelihoodSelect(BINOMIAL);
		
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
		
		
//		aMatrix = new Haplotypes(aMap, 5);
//		for (int i = 0; i < aMatrix.getHaplotypesCount(); i++) {
//			aMatrix.randomSeq(i);
//			System.out.println(aMatrix.getHaplotype(i) );
//		}

		srL = new ShortReadLikelihood(aMap, aMatrix);
		System.out.println("Likelihood:\t"+srL.getLogLikelihood()+"\n");
		
		double likelihood = srL.getLogLikelihood(); 
		
		MCMCCriterion criterion = new MCMCCriterion();
        double hastingsRatio = 0.0;
        double[] logr = {-Double.MAX_VALUE};
        
		for (int i = 0; i < 1e4; i++) {
			
			aMatrix.swapBase();
			srL.updateHaplotypes(aMatrix);
//			System.out.print(srL.getLogLikelihood() +"\t" );
//			srL = new ShortReadLikelihood(aMap, aMatrix);
//			System.out.println(srL.getLogLikelihood());
			double newL = srL.getLogLikelihood();
			
			
			
			boolean accept = criterion.accept(likelihood, newL, hastingsRatio, logr);
//			System.out.println(likelihood +"\t"+newL +"\t"+ logr[0] +"\t"+  accept);
			if(accept){
				likelihood = newL;
				srL.acceptState();
			}
			else{
				aMatrix.reject();
				srL.restoreState();
				
			}
			
				
		}
		for (int i = 0; i < aMatrix.getHaplotypesCount(); i++) {
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
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		seqs = new String[]{"TTGGCCAA"};
		SimpleAlignment hap = AlignmentUtils.createAlignment(seqs);
		
		Haplotypes aMatrix = new Haplotypes(aMap, hap);
		ShortReadLikelihood srL = new ShortReadLikelihood(aMap, aMatrix);
		
		double logLikelihood = srL.calculateLogLikelihoodSelect(BINOMIAL);
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

		String[] seqs = new String[]{
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		String[] haps = new String[]{
				"AAAACCCCGGGGTTTT",
				"ACGTACGTACGTACGT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		Haplotypes aMatrix = new Haplotypes(aMap, hapAlignment);

		ShortReadLikelihood srL = new ShortReadLikelihood(aMap, aMatrix);
		double logLikelihood = srL.getLogLikelihood();

		
		double expected = -0.0819407854607265; // log(sum(dpois(seq(0,8,by=2),8*E)))
		

		assertEquals(expected, logLikelihood, 1e-10);
		logLikelihood = srL.calculateLogLikelihoodSelect(BINOMIAL);
		expected = -7.6461084423155947821e-09; //  log(sum(dbinom(0:4,8,E)))
		assertEquals(expected, logLikelihood, 1e-10);

		
	}
	
	@Test
	public void testRandomQuickRun(){

		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/";

		String trueAlignmentFile = "1110_10_org_3.phyml";
		String shortReadFile = "1110_10_align_100.fasta";
		
		
		DataImporter dataImporter = new DataImporter(dataDir);
		
		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
		
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		Haplotypes trueHaplotypes = new Haplotypes(aMap, trueAlignment);
		
		ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(aMap, trueHaplotypes);

		
		System.out.println("Likelihood:\t"+srpLikelihood.getLogLikelihood() 
				+"\t"+ trueHaplotypes.calculateSPS());
		int[][] sps = Haplotypes.calculeteSPSArray(trueHaplotypes, trueHaplotypes);
		for (int i = 0; i < sps.length; i++) {
				System.out.println(Arrays.toString(sps[i]));
		}
		

		Haplotypes haplotypes = new Haplotypes(aMap, trueAlignment.getSequenceCount());
		for (int i = 0; i < haplotypes.getHaplotypesCount(); i++) {
			System.out.println(haplotypes.getHaplotype(i) );
			System.out.println(trueHaplotypes.getHaplotype(i) );
			System.out.println();
		}

		srpLikelihood = new ShortReadLikelihood(aMap, haplotypes);
		
		System.out.println("Likelihood:\t"+srpLikelihood.getLogLikelihood() 
				+"\t"+ haplotypes.calculateSPS());
		sps = Haplotypes.calculeteSPSArray(trueHaplotypes, haplotypes);
		for (int i = 0; i < sps.length; i++) {
				System.out.println(Arrays.toString(sps[i]));
		}
		
		double likelihood = srpLikelihood.getLogLikelihood(); 
	
        
		MCMCCriterion criterion = new MCMCCriterion();
        double hastingsRatio = 0.0;
        double[] logr = {-Double.MAX_VALUE};
        
        
        
        int thinning = 1000;
        for (int i = 0; i < 1e2; i++) {

			haplotypes.swapBase();
			srpLikelihood.updateHaplotypes(haplotypes);

			double newL = srpLikelihood.getLogLikelihood();

			boolean accept = criterion.accept(likelihood, newL, hastingsRatio, logr);
//				System.out.println(likelihood +"\t"+newL +"\t"+ logr[0] +"\t"+  accept);
			if(accept){
				likelihood = newL;
				srpLikelihood.acceptState();
			}
			else{
				haplotypes.reject();
				srpLikelihood.restoreState();
			}
			
					
		}

		System.out.println("===== Done =====");
		System.out.println("Likelihood:\t"+srpLikelihood.getLogLikelihood() 
				+"\t"+ haplotypes.calculateSPS());
		sps = Haplotypes.calculeteSPSArray(trueHaplotypes, haplotypes);
		for (int i = 0; i < sps.length; i++) {
				System.out.println(Arrays.toString(sps[i]));
		}
		

	
	}
//
//	@Test
//	public void testCalculateLikelihood3() {
//		
//		
//		shortReads.add("AAAACCCC");
//		haplotypes.add("AAAACCCCAAAA"); //0 2 4 6 8 mismatches
//		haplotypes.add("AAAACCCCCCCC"); // 0 1 2 3 4 mismatches
//		ShortReadLikelihood srL = new ShortReadLikelihood(shortReads, haplotypes);
//		
//		double logLikelihood = srL.calculateLogLikelihoodSelect(BINOMIAL);
//		double expected = 0.653015821111653; // log(sum(dpois(c(seq(0,8,by=2),0:4),8*E)))
//		assertEquals("multiple haplotypes", expected, logLikelihood, 1e-10);
//		
//		logLikelihood = srL.calculateLogLikelihoodSelect(BINOMIAL);
//		expected =  0.65260860360614092457; // log(sum(dbinom(c(seq(0,8,by=2),0:4),8,E)))
//		assertEquals("multiple haplotypes", expected, logLikelihood, 1e-10);
//		
//		
//		shortReads.add("AACCCAA");	//		4 3 1 1 3 4 mismatches
//		 							// and	4 3 2 3 4 4 mismatches
//		srL.setShortRead(shortReads);
//		logLikelihood = srL.calculateLogLikelihoodSelect(BINOMIAL);
//		expected  =  -1.299909707366587952; 
//		// log(sum(dpois(c(seq(0,8,by=2),0:4),8*E))) + log(sum(dpois(c(4,3,1,1,3,4,4,3,2,3,4,4),7*E)))
//		assertEquals("multiple short reads", expected, logLikelihood, 1e-10);
//		
//		logLikelihood = srL.calculateLogLikelihoodSelect(BINOMIAL);
//		expected =  -1.2931321407960054692;
//		// log(sum(dbinom(c(seq(0,8,by=2),0:4),8,E))) + log(sum(dbinom(c(4,3,1,1,3,4,4,3,2,3,4,4),7,E)))
//		assertEquals("multiple haplotypes", expected, logLikelihood, 1e-10);
//		
//		haplotypes.add("TTTTGGGGATTGGACAC");
////		                          8,8,8,8,8,7,7,6,6,6	mismatches
////		                          7,7,6,6,7,7,7,5,5,6,5 mismatches
//		srL.setHaplotypes(haplotypes);
//		logLikelihood = srL.calculateLogLikelihoodSelect(BINOMIAL);
//		expected  = -1.2999093162319743655; //log(sum(dpois(c(seq(0,8,by=2),0:4, 8,8,8,8,8,7,7,6,6,6),8*E)))+ log(sum(dpois(c(4,3,1,1,3,4,4,3,2,3,4,4,  7,7,6,6,7,7,7,5,5,6,5),7*E)))
//		assertEquals(expected, logLikelihood, 1e-10);
//
//		logLikelihood = srL.calculateLogLikelihoodSelect(BINOMIAL);
//		expected =  -1.2931320799871470761;
//		// log(sum(dpois(c(seq(0,8,by=2),0:4, 8,8,8,8,8,7,7,6,6,6),8*E)))+ log(sum(dpois(c(4,3,1,1,3,4,4,3,2,3,4,4,  7,7,6,6,7,7,7,5,5,6,5),7*E)))
//		assertEquals("multiple haplotypes", expected, logLikelihood, 1e-10);
//
//		
//	}
//	
//	@Test
//	public void testCompoundLikelihood() {
//
//		
//
//		shortReads.add("AAAACCCC");
//		shortReads.add("AACCCAA");
//		
//		haplotypes.add("AAAACCCCAAAA"); //0 2 4 6 8 mismatches
//		haplotypes.add("AAAACCCCCCCC"); // 0 1 2 3 4 mismatches
//		haplotypes.add("TTTTGGGGATTGGACAC");
//		
//		ShortReadLikelihood srL = new ShortReadLikelihood(shortReads, haplotypes);
//		
//		double logLikelihood = srL.getLogLikelihood();
//		double expected =  -1.2999093162319743655; //log(sum(dpois(c(seq(0,8,by=2),0:4, 8,8,8,8,8,7,7,6,6,6),8*E)))+ log(sum(dpois(c(4,3,1,1,3,4,4,3,2,3,4,4,  7,7,6,6,7,7,7,5,5,6,5),7*E)))
//		assertEquals("multiple short reads", expected, logLikelihood, 1e-10);
//		
//        List<Likelihood> likelihoods = new ArrayList<Likelihood>();        
//        likelihoods.add(srL);
//        Likelihood CSRLikelihood = new CompoundLikelihood(0, likelihoods);
//        CSRLikelihood.setId(ShortReadLikelihood.NAME);
////        System.out.println(((CompoundLikelihood) CSRLikelihood).getDiagnosis());
////        System.out.println(((CompoundLikelihood) CSRLikelihood).getReport());
////        System.out.println(CSRLikelihood.getId());
////        System.out.println(CSRLikelihood.prettyName());
////        System.out.println(CSRLikelihood.getModel().getModelName());
//        logLikelihood = CSRLikelihood.getLogLikelihood();
//        assertEquals("multiple short reads", expected, logLikelihood, 1e-10);
//        
//
//	}
//	
//	@Test
//	public void testDifferentImplementation100Reads(){
//		
//		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/Stage0/";
//		
//		if( Files.exists(Paths.get(dataDir)) ){
//	
//			String trueAlignmentFile = "121101_true_seqs.fasta";
//			String shortReadFile = "121101_short_reads_100.fasta";
//			
//			DataImporter dataImporter = new DataImporter(dataDir);
//			Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
//			Sequences shortReads = dataImporter.importSequence(shortReadFile);
//			
//			ShortReadLikelihood srl = new ShortReadLikelihood(shortReads, trueAlignment);
//		
//			double expected = -106.642313524702; //100 reads
//			
//			// calculateShoreReadLikelihood1 ~ 14500-15000
//			// calculateShoreReadLikelihood2 ~ 800-900
//			
//			double actual = srl.calculateLogLikelihoodSelect(0);
//			assertEquals("Compare to python version, 100 reads", expected, actual, 1e-10);
//			actual = srl.calculateLogLikelihoodSelect(1);
//			assertEquals("Compare to python version, 100 reads", expected, actual, 1e-10);
////			actual = srl.getLogLikelihoodSelect(2);
////			assertEquals("Compare to python version, 100 reads", expected, actual, 1e-10);
//		
//			
//		} else {
//			fail("change hardcoded dataDir: " + dataDir);
//		}
//	}
//
//	
//
//	@Test
//	public void testDifferentImplementation10Reads(){
//		
//		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/Stage0/";
//		
//		if( Files.exists(Paths.get(dataDir)) ){
//	
//			String trueAlignmentFile = "121101_true_seqs.fasta";
//			String shortReadFile = "121101_short_reads_10.fasta";
//			
//			DataImporter dataImporter = new DataImporter(dataDir);
//			Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
//			Sequences shortReads = dataImporter.importSequence(shortReadFile);
//			
//			ShortReadLikelihood srl = new ShortReadLikelihood(shortReads, trueAlignment);
//		
//			double expected = -11.85523370076271; //10  reads
//			
//			double actual = srl.calculateLogLikelihoodSelect(0);
//			assertEquals("Compare to python version, 10 reads", expected, actual, 1e-10);
//			actual = srl.calculateLogLikelihoodSelect(1);
//			assertEquals("Compare to python version, 10 reads", expected, actual, 1e-10);
////			actual = srl.getLogLikelihoodSelect(2);
////			assertEquals("Compare to python version, 100 reads", expected, actual, 1e-10);
//		
//			
//		} else {
//			fail("change hardcoded dataDir: " + dataDir);
//		}
//	}

}
