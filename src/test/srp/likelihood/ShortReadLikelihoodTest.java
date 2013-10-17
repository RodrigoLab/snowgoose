package test.srp.likelihood;

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
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.HaplotypeModelUtils;
import srp.haplotypes.Operation;
import srp.haplotypes.operator.AbstractMultiBasesOperator;
import srp.haplotypes.operator.AbstractSingleBaseOperator;
import srp.haplotypes.operator.ColumnOperator;
import srp.haplotypes.operator.HaplotypeRecombinationOperator;
import srp.haplotypes.operator.HaplotypeSwapSectionOperator;
import srp.haplotypes.operator.MultiBasesEmpiricalOperator;
import srp.haplotypes.operator.MultiBasesOperator;
import srp.haplotypes.operator.MultiBasesUniformOperator;
import srp.haplotypes.operator.SingleBaseEmpiricalOperator;
import srp.haplotypes.operator.SingleBaseFrequencyOperator;
import srp.haplotypes.operator.SingleBaseOperator;
import srp.haplotypes.operator.SingleBaseUniformOperator;
import srp.likelihood.ShortReadLikelihood;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCCriterion;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.math.MathUtils;

public class ShortReadLikelihoodTest {

	
	private HaplotypeModel haplotypeModelH4;
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {

	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {


		String dataDir = "/home/sw167/workspaceSrp/ABI/unittest/";
	
		String trueAlignmentFile = "H4_haplotypes.phyml";
		String shortReadFile = "H4_srp.fasta";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		
		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
		
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		haplotypeModelH4 = new HaplotypeModel(aMap, trueAlignment.getSequenceCount());
	
	}

	@After
	public void tearDown() throws Exception {
	}

	
	
	@Test
	public void testCalculateLikelihoodIdentical() {

		String[] seqs = new String[]{"AACCGGTT"};
	
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
//		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		
		HaplotypeModel haplotypeModel = AlignmentUtils.createHaplotypeModel(seqs);
		ShortReadLikelihood srL = new ShortReadLikelihood(haplotypeModel);
		
		double logLikelihood = srL.getLogLikelihood();
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
		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		
		haplotypeModel = new HaplotypeModel(aMap, alignment);
		srL = new ShortReadLikelihood(haplotypeModel);
		
		
		logLikelihood = srL.getLogLikelihood();
		expected = -0.086061253223681313806*4; //dbinom(0,8,E,log=T)
		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);
		
	}



	@Test
	public void testCalculateLikelihoodSingle() throws Exception {
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
//				"AAAAACCCCCGGGGGTTTTTACGTACACTATATATA",
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//				"AAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, hapAlignment);

		ShortReadLikelihood srL = new ShortReadLikelihood(haplotypeModel);
		List<Likelihood> likelihoods = new ArrayList<Likelihood>();

		
		likelihoods.add(srL);
		Likelihood shortReadlikelihood = new CompoundLikelihood(-1, likelihoods);
		
        SingleBaseOperator op = new SingleBaseOperator(haplotypeModel, 0);
        OperatorSchedule schedule = new SimpleOperatorSchedule();
		schedule.addOperator(op);
        
		int length = 1000;
		MCMCOptions options = new MCMCOptions();
		options.setChainLength(length);
//		options.setUseCoercion(false); // autoOptimize = true
		options.setFullEvaluationCount(length);

		MCMC mcmc = new MCMC("mcmc1");
		mcmc.init(options, shortReadlikelihood, schedule, null);
		mcmc.run();
		
	}
	
	@Test
	public void testCalculateLikelihoodSomething() throws Exception {
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
//				"AAAAACCCCCGGGGGTTTTTACGTACACTATATATA",
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//				"AAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, hapAlignment);

		ShortReadLikelihood likelihood = new ShortReadLikelihood(haplotypeModel);

		// from MarkovChain.runChain();		
        likelihood.makeDirty();
        double currentScore = likelihood.getLogLikelihood();
        long currentState = 0;
        Model currentModel = likelihood.getModel(); //should be posterior
        MCMCCriterion acceptor = new MCMCCriterion();
        MCMCOperator op = new SingleBaseOperator(haplotypeModel, 0);
        
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
				"AAAAAAAAATGTGTTTT....",
				".....CCCCCCCCCCCCCCCCCCCTTTTCCCC....",
				"..........GGGGGGGGGGGGGGCGCGTATAGGGG",
				"...............TTTTTTTTTACACTATA....",
//    				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//    				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG"
//    				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		String[] haps = new String[]{
//    				"AAAAACCCCCGGGGGTTTTTACGTACACTATATATA",
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//    				"AAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, hapAlignment);

		ShortReadLikelihood srL = new ShortReadLikelihood(haplotypeModel);

		MCMCCriterion acceptor = new MCMCCriterion();
		double[] logr = { 0.0 };

		SingleBaseOperator op = new SingleBaseOperator(haplotypeModel, 0);
		MCMCOperator mcmcOperator = new SingleBaseOperator(haplotypeModel, 0);

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
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		seqs = new String[]{"TTGGCCAA"};
		SimpleAlignment hap = AlignmentUtils.createAlignment(seqs);
		
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, hap);
		ShortReadLikelihood srL = new ShortReadLikelihood(haplotypeModel);
		
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
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(srps);
		
		String[] haps = new String[]{
				"AAAACCCCGGGGTTTT",
				"ACGTACGTACGTACGT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, hapAlignment);

		ShortReadLikelihood srL = new ShortReadLikelihood(haplotypeModel);
		double logLikelihood = srL.getLogLikelihood();

		
		double expected = -6.4817467758693254609;
//		E = 0.0107; log(sum(dbinom(c(0,6),8,E))) + log(sum(dbinom(c(3,2),6,E)))

		assertEquals(expected, logLikelihood, 1e-10);

		
	}
	
	@Test
	public void testRandomQuickRun() throws Exception{

		String dataDir = "/home/sw167/workspaceSrp/ABI/unittest/";

		String trueAlignmentFile = "H4_haplotypes.phyml";
		String shortReadFile = "H4_srp.fasta";
		
		
		DataImporter dataImporter = new DataImporter(dataDir);
		
		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
		
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
//		HaplotypeModel trueHaplotypes = new HaplotypeModel(aMap, trueAlignment);
		
//		ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(aMap, trueHaplotypes);

		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, trueAlignment.getSequenceCount());

		ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(haplotypeModel);
		
		double likelihood = srpLikelihood.getLogLikelihood(); 
		double initLikelihood = likelihood;
        
		MCMCCriterion criterion = new MCMCCriterion();
        double hastingsRatio = 0.0;
        double[] logr = {-Double.MAX_VALUE};
        
        SingleBaseOperator op = new SingleBaseOperator(haplotypeModel, 0);
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
				haplotypeModel.reject();
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
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(srps);
		
		String[] haps = new String[]{
				"AAAACCCCGGGGTTTT",
				"ACGTACGTACGTACGT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, hapAlignment);

		ShortReadLikelihood srL = new ShortReadLikelihood(haplotypeModel);
		
		double logLikelihood = srL.getLogLikelihood();

		
		double expected = -6.4817467758693254609;
//		E = 0.0107; log(sum(dbinom(c(0,6),8,E))) + log(sum(dbinom(c(3,2),6,E)))

		assertEquals(expected, logLikelihood, 1e-10);

        List<Likelihood> likelihoods = new ArrayList<Likelihood>();        
        likelihoods.add(srL);
        
        Likelihood CSRLikelihood = new CompoundLikelihood(0, likelihoods);
        CSRLikelihood.setId(ShortReadLikelihood.NAME);

        logLikelihood = CSRLikelihood.getLogLikelihood();
        assertEquals(expected, logLikelihood, 1e-10);
        

	}

	

	@Test
	public void testCalculateLikelihoodTime() throws Exception {

//		String[] seqs = new String[]{
//				"AAAAAATTTTT.........",
//				"..CCCCCCCCCCCCC.....",
//				"......GGGGGGGGGGG...",
//				"........TTTTTTTTT...",
//				"............ACGTACGT",
//				};
//		
//		String[] haps = new String[]{
//
//				"CCCCCTTTTTAAAAAGGGGG",
//				"ACGTACGTACGTACGTACGT",
//
//				};
//		
//		HaplotypeModel haplotypeModel = AlignmentUtils.createHaplotypeModel(seqs, haps);
//
//		ShortReadLikelihood srL = new ShortReadLikelihood(haplotypeModel);
//
		
		

		String dataDir = "/home/sw167/workspaceSrp/ABI/unittest/";

		String trueAlignmentFile = "H4_haplotypes.phyml";
		String shortReadFile = "H4_srp.fasta";
		
		
		DataImporter dataImporter = new DataImporter(dataDir);
		
		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
		
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);

		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, trueAlignment.getSequenceCount());

		ShortReadLikelihood srL = new ShortReadLikelihood(haplotypeModel);
		
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
			
			haplotypeModel.swapHaplotypeColumn(new int[]{pos, 'T'});
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
		
	@Test
	public void testCalculateLikelihoodColumn() throws Exception {
		
		Parameter freqs = new Parameter.Default("frequency", haplotypeModelH4.getStateFrequencies());
		MCMCOperator op = new ColumnOperator(haplotypeModelH4, haplotypeModelH4.getHaplotypeCount(), freqs, null);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, ColumnOperator.OP);
	}
	@Test
	public void testCalculateLikelihoodSingleBase() throws Exception {
		
		MCMCOperator op = new SingleBaseOperator(haplotypeModelH4, 0);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, AbstractSingleBaseOperator.OP);
		
	}
	@Test
	public void testCalculateLikelihoodSingleBaseUniform() throws Exception {
		
		MCMCOperator op = new SingleBaseUniformOperator(haplotypeModelH4, 0);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, AbstractSingleBaseOperator.OP);
	}
	@Test
	public void testCalculateLikelihoodSingleBaseEmpirical() throws Exception {
		
		MCMCOperator op = new SingleBaseEmpiricalOperator(haplotypeModelH4, 0);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, AbstractSingleBaseOperator.OP);
	}
	
	@Test
	public void testCalculateLikelihoodSingleBaseFrequency() throws Exception {
			

		Parameter freqs = new Parameter.Default("frequency", haplotypeModelH4.getStateFrequencies());
		MCMCOperator op = new SingleBaseFrequencyOperator(haplotypeModelH4, freqs);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, AbstractSingleBaseOperator.OP);
		
		
	}
	@Test
	public void testCalculateLikelihoodMultiBases() throws Exception {
		
		MCMCOperator op = new MultiBasesOperator(haplotypeModelH4, 50, CoercionMode.COERCION_OFF);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, AbstractMultiBasesOperator.OP);
	}
	@Test
	public void testCalculateLikelihoodMultiBasesEmpirical() throws Exception {
		
		MCMCOperator op = new MultiBasesEmpiricalOperator(haplotypeModelH4, 50, CoercionMode.COERCION_OFF);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, AbstractMultiBasesOperator.OP);
	}
	@Test
	public void testCalculateLikelihoodMultiBasesUniform() throws Exception {
		
		MCMCOperator op = new MultiBasesUniformOperator(haplotypeModelH4, 50, CoercionMode.COERCION_OFF);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, AbstractMultiBasesOperator.OP);
	}
	@Test
	public void testCalculateLikelihoodSwapSectionRecombination() throws Exception {
		
		MCMCOperator op = new HaplotypeRecombinationOperator(haplotypeModelH4, 0);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, HaplotypeRecombinationOperator.OP);
	}
	@Test
	public void testCalculateLikelihoodSwapSection() throws Exception {
		
		MCMCOperator op = new HaplotypeSwapSectionOperator(haplotypeModelH4, 50, null);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, HaplotypeSwapSectionOperator.OP);
	}
	
	public void runTestCalculateSrpLikelihoodOperators(HaplotypeModel haplotypeModel, MCMCOperator op, Object expectedOperation) throws Exception {
	
		
		ShortReadLikelihood srL = new ShortReadLikelihood(haplotypeModel);
		assertEquals(Operation.NONE, srL.getOperation());
		assertEquals(Operation.NONE, haplotypeModel.getSwapInfo().getOperation());
		for (int i = 0; i < 100; i++) {
	    	srL.storeModelState();
	        op.operate();
	
	        double score = srL.getLogLikelihood();
	        assertEquals(expectedOperation, srL.getOperation());
	        assertEquals(expectedOperation, haplotypeModel.getSwapInfo().getOperation());
	        
	        HaplotypeModel duplicateHaplotypeModel = HaplotypeModelUtils.copyHaplotypeModel(haplotypeModel);
	        ShortReadLikelihood srLFull = new ShortReadLikelihood(duplicateHaplotypeModel);
	        double expected = srLFull.getLogLikelihood();
	        
	        assertEquals(Operation.NONE, srLFull.getOperation());
	        assertEquals(Operation.NONE, duplicateHaplotypeModel.getSwapInfo().getOperation());
	        assertEquals(expected, score, 0);
	
	        boolean accept = MathUtils.nextBoolean();    		
	        if (accept) {
	            op.accept(0);
	            srL.acceptModelState();
	        } else {
	            op.reject();
	            srL.restoreModelState();
	        }
	    }
		
		
	}
}
