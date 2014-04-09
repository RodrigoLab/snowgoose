package test.srp.likelihood.haplotypes;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;

import srp.evolution.shortreads.ShortReadMapping;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.HaplotypeModel;
import srp.likelihood.haplotypes.ShortReadsHaplotypeLikelihood;
import srp.operator.haplotypes.BaseSingleOperator;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCCriterion;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.math.MathUtils;

public class ShortReadsHaplotypeLikelihoodTest {

	
	private HaplotypeModel haplotypeModel;
	private ShortReadMapping  srpMap;
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {

	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {


		String dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/";
	
		String trueAlignmentFile = "H4_haplotypes.phyml";
		String shortReadFile = "H4_srp.fasta";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		
		Alignment shortReads = dataImporter.importShortReads(shortReadFile);
		srpMap = new ShortReadMapping(shortReads);
		
		
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		haplotypeModel = new HaplotypeModel(trueAlignment);
	
	}

	@After
	public void tearDown() throws Exception {
	}

	
	
	@Test
	public void testCalculateLikelihoodIdentical() {

		String[] seqs = new String[]{"AACCGGTT"};
	
		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
		
		
//		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
//		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
//		ShortReadMapping srpMap = ShortReadMapping
		HaplotypeModel haplotypeModel = new HaplotypeModel(alignment);
		ShortReadsHaplotypeLikelihood srL = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		
		double logLikelihood = srL.getLogLikelihood();
		double expected = -0.086061253223681313806; //dbinom(0,8,E,log=T)
		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);
		
		seqs = new String[]{
				".....AACCGGTT........",
				".....AACCGGTT........",
				".....AACCGGTT........",
				".....AACCGGTT........",

				};
		alignment = AlignmentUtils.createAlignment(seqs);
		srpMap = new ShortReadMapping(alignment);
		seqs = new String[]{"WWWWWAACCGGTTYYYYYYYY"};
		
		alignment = AlignmentUtils.createAlignment(seqs);
		
		haplotypeModel = new HaplotypeModel(alignment);
		srL = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		
		
		logLikelihood = srL.getLogLikelihood();
		expected = -0.086061253223681313806*4; //dbinom(0,8,E,log=T)
		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);
		
	}



	@Test
	public void testCalculateLikelihoodSingle() throws Exception {
		String[] seqs = new String[]{
				"AAAAAAAAATGTGTTTT...................",
				".....CCCCCCCCCCCCCCCCCCCTTTTCCCC....",
				"..........GGGGGGGGGGGGGGCGCGTATAGGGG",
				"...............TTTTTTTTTACACTATA....",
//				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG"
//				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG
				};
		Alignment alignment = AlignmentUtils.createAlignment(seqs);
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
		
		String[] haps = new String[]{
//				"AAAAACCCCCGGGGGTTTTTACGTACACTATATATA",
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//				"AAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(hapAlignment);
		ShortReadsHaplotypeLikelihood srL = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		
		List<Likelihood> likelihoods = new ArrayList<Likelihood>();

		
		likelihoods.add(srL);
		Likelihood shortReadlikelihood = new CompoundLikelihood(-1, likelihoods);
		
        BaseSingleOperator op = new BaseSingleOperator(haplotypeModel);
        OperatorSchedule schedule = new SimpleOperatorSchedule();
		schedule.addOperator(op);
        
		int length = 1000;
		MCMCOptions options = new MCMCOptions(length);
//		options.setChainLength(length);
//		options.setUseCoercion(false); // autoOptimize = true
//		options.setFullEvaluationCount(length);

		MCMC mcmc = new MCMC("mcmc1");
		mcmc.init(options, shortReadlikelihood, schedule, null);
		mcmc.run();
		
	}
	
	@Test
	public void testCalculateLikelihoodSomething() throws Exception {
		String[] seqs = new String[]{
				"AAAAAAAAATGTGTTTT...................",
				".....CCCCCCCCCCCCCCCCCCCTTTTCCCC....",
				"..........GGGGGGGGGGGGGGCGCGTATAGGGG",
				"...............TTTTTTTTTACACTATA....",
//				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG"
//				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG
				};
		
		Alignment alignment = AlignmentUtils.createAlignment(seqs);
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
		
		String[] haps = new String[]{
//				"AAAAACCCCCGGGGGTTTTTACGTACACTATATATA",
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//				"AAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(hapAlignment);
		ShortReadsHaplotypeLikelihood likelihood = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		
		

		// from MarkovChain.runChain();		
        likelihood.makeDirty();
        double currentScore = likelihood.getLogLikelihood();
        long currentState = 0;
        Model currentModel = likelihood.getModel(); //should be posterior
        MCMCCriterion acceptor = new MCMCCriterion();
        MCMCOperator op = new BaseSingleOperator(haplotypeModel);
        
        double[] logr = {0.0};
        while ( currentState < 100) {

            final MCMCOperator mcmcOperator = op;
            double oldScore = currentScore;
            currentModel.storeModelState();
            
            double hastingsRatio = 1.0;
            logr[0] = -Double.MAX_VALUE;
            hastingsRatio = mcmcOperator.operate();

            double score = likelihood.getLogLikelihood();
            boolean accept = acceptor.accept(oldScore, score, hastingsRatio, logr);
            double deviation = score - oldScore;

            if (accept) {
                mcmcOperator.accept(deviation);
                currentModel.acceptModelState();
                currentScore = score;

            } else {
                mcmcOperator.reject();
                currentModel.restoreModelState();
            }
            currentState += 1;
        }
	}

	@Test
	public void testCalculateLikelihoodSomething2() throws Exception {
		String[] seqs = new String[]{
				"AAAAAAAAATGTGTTTT...................",
				".....CCCCCCCCCCCCCCCCCCCTTTTCCCC....",
				"..........GGGGGGGGGGGGGGCGCGTATAGGGG",
				"...............TTTTTTTTTACACTATA....",
//    				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//    				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG"
//    				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG
				};
		Alignment alignment = AlignmentUtils.createAlignment(seqs);
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
		
		String[] haps = new String[]{
//				"AAAAACCCCCGGGGGTTTTTACGTACACTATATATA",
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//				"AAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(hapAlignment);
		ShortReadsHaplotypeLikelihood srL = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		
		MCMCCriterion acceptor = new MCMCCriterion();
		double[] logr = { 0.0 };

		BaseSingleOperator op = new BaseSingleOperator(haplotypeModel);
		MCMCOperator mcmcOperator = new BaseSingleOperator(haplotypeModel);

		double currentScore = srL.getLogLikelihood();
		
		for (int i = 0; i < 1e4; i++) {
			
			srL.storeModelState();
			double oldScore = currentScore;
			logr[0] = -Double.MAX_VALUE;

			double hastingsRatio = mcmcOperator.operate();
//			double hastingsRatio = op.doOperation();

			double score = srL.getLogLikelihood();
			boolean accept = acceptor.accept(oldScore, score, hastingsRatio, logr);
            double deviation = score - oldScore;
            
			if (accept) {
				mcmcOperator.accept(deviation);
				srL.acceptModelState();
				currentScore = score;
			} else {
				mcmcOperator.reject();
				srL.restoreModelState();
			}
		}

	}
	@Test
	public void testCalculateLikelihoodAllDiff() {

		String[] seqs = new String[]{"AACCGGTT"};
		Alignment alignment = AlignmentUtils.createAlignment(seqs);
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
		
		
		String[] haps = new String[]{"TTGGCCAA"};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(hapAlignment);
		ShortReadsHaplotypeLikelihood srL = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);

		double logLikelihood = srL.getLogLikelihood();
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

		String[] srps = new String[]{
//				"AAAACCCCGGGGTTTT",
//				"ACGTACGTACGTACGT",
				"....CCCCGGGG....", // 0/8 + 6/8
				"..AAACGT........", // 3/6 + 2/6
				};
		Alignment alignment = AlignmentUtils.createAlignment(srps);
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
		
		
		String[] haps = new String[]{
				"AAAACCCCGGGGTTTT",
				"ACGTACGTACGTACGT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(hapAlignment);
		ShortReadsHaplotypeLikelihood srL = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);

		double logLikelihood = srL.getLogLikelihood();

		
		double expected = -6.4817467758693254609;
//		E = 0.0107; log(sum(dbinom(c(0,6),8,E))) + log(sum(dbinom(c(3,2),6,E)))

		assertEquals(expected, logLikelihood, 1e-10);

		
	}
	
	@Test
	public void testRandomQuickRun() throws Exception{

		ShortReadsHaplotypeLikelihood srpLikelihood = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		
		double likelihood = srpLikelihood.getLogLikelihood(); 
		double initLikelihood = likelihood;
        
		MCMCCriterion criterion = new MCMCCriterion();
        double hastingsRatio = 0.0;
        double[] logr = {-Double.MAX_VALUE};
        
        BaseSingleOperator op = new BaseSingleOperator(haplotypeModel);
        for (int i = 0; i < 1e2; i++) {
        	op.doOperation();
        	srpLikelihood.storeModelState();
			srpLikelihood.makeDirty();

			double newL = srpLikelihood.getLogLikelihood();

			boolean accept = criterion.accept(likelihood, newL, hastingsRatio, logr);
			if(accept){
				likelihood = newL;
//				srpLikelihood.acceptState();
			}
			else{
//				haplotypeModel.resreject();
				srpLikelihood.restoreModelState();
			}
			
			
		}
        double finalLogLikelihood = srpLikelihood.getLogLikelihood(); 
		assertTrue(finalLogLikelihood >= initLikelihood);

		

	
	}
	
	@Test
	public void testCompoundLikelihood() {


		String[] srps = new String[]{
//				"AAAACCCCGGGGTTTT",
//				"ACGTACGTACGTACGT",
				"....CCCCGGGG....", // 0/8 + 6/8
				"..AAACGT........", // 3/6 + 2/6
				};
		Alignment alignment = AlignmentUtils.createAlignment(srps);
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
		
		
		String[] haps = new String[]{
				"AAAACCCCGGGGTTTT",
				"ACGTACGTACGTACGT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(hapAlignment);
		ShortReadsHaplotypeLikelihood srL = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		
		double logLikelihood = srL.getLogLikelihood();

		
		double expected = -6.4817467758693254609;
//		E = 0.0107; log(sum(dbinom(c(0,6),8,E))) + log(sum(dbinom(c(3,2),6,E)))

		assertEquals(expected, logLikelihood, 1e-10);

        List<Likelihood> likelihoods = new ArrayList<Likelihood>();        
        likelihoods.add(srL);
        
        Likelihood CSRLikelihood = new CompoundLikelihood(0, likelihoods);
//        CSRLikelihood.setId(ShortReadLikelihood.NAME);

        logLikelihood = CSRLikelihood.getLogLikelihood();
        assertEquals(expected, logLikelihood, 1e-10);
        

	}

	

	@Test
	public void testCalculateLikelihoodTime() throws Exception {

		ShortReadsHaplotypeLikelihood srL = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		
		double logLikelihood = srL.getLogLikelihood();
		long time1 = System.currentTimeMillis();
		double ite =  1e3;
		for (int i = 0; i < ite; i++) {
			srL.makeDirty();
//			logLikelihood = srL.getLogLikelihood();
		}
		long time2 = System.currentTimeMillis();
		System.out.println((time2 - time1)/ite + " ms per full calculation");
		
		
		ite = 1e5;
		time1 = System.currentTimeMillis();
		long total = 0;
		for (int i = 0; i < ite; i++) {

			srL.storeModelState();
//			haplotypeModel.swapHaplotypeSingleBase(Operation.SWAPSINGLE, new int[]{600,'T'});
			int pos = MathUtils.nextInt(haplotypeModel.getHaplotypeLength());
			
//			haplotypeModel.swapHaplotypeColumn(new int[]{pos, 'T'});
			srL.makeDirty();

			time1 = System.currentTimeMillis();
			logLikelihood = srL.getLogLikelihood();
			
			
			time2 = System.currentTimeMillis();
			total += (time2-time1);
			srL.restoreModelState();
		}
		time2 = System.currentTimeMillis();
		
		System.out.println(total/ite + " ms per calculations");
		
		
	}

}
