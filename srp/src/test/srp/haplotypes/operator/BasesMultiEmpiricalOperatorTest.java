package test.srp.haplotypes.operator;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.likelihood.ShortReadLikelihood;
import srp.haplotypes.operator.BasesMultiEmpiricalOperator;
import dr.evolution.alignment.SimpleAlignment;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.trace.Trace;

public class BasesMultiEmpiricalOperatorTest {

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
	public void testGetOperatorName() {
		String[] seqs = new String[]{
				"AAAAAAAAATGTGTTTT....",
				".....CCCCCCCCCCCCCCCCCCCTTTTCCCC....",
				"..........GGGGGGGGGGGGGGCGCGTATAGGGG",
				"...............TTTTTTTTTACACTATA....",
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, 3);

		int nBases = 10;
		CoercableMCMCOperator operator = new BasesMultiEmpiricalOperator(haplotypeModel, nBases, CoercionMode.COERCION_OFF);
    	assertEquals(operator.getOperatorName(), "BasesMultiEmpiricalOperator");
    	assertEquals(operator.getRawParameter(), nBases, 0);
    	assertEquals(operator.getCoercableParameter(), Math.log(nBases-1), 1e-10); 
    	
    	assertEquals(operator.getPerformanceSuggestion(), "");
	}

	
	@Test
	public void testDoOperation() throws OperatorFailedException {
		String[] seqs = new String[]{
				"GGGGGGGGGGGGG.....",
				".....CCCCCCCCCCCCCCCCCCCCCCCCCCC....",
				"..........GGGGGGGGGGGGGGCGCGGGGGGGGG",
//				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG"
//				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG
				};
		
		String[] haps = new String[]{
//				"AAAAACCCCCGGGGGTTTTTACGTACACTATATATA"
//				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
				"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
				};
		
		HaplotypeModel haplotypeModel = AlignmentUtils.createHaplotypeModel(seqs, haps);
    	BasesMultiEmpiricalOperator operator = new BasesMultiEmpiricalOperator(haplotypeModel, 5, CoercionMode.COERCION_OFF);
    	
    	
    	for (int i = 0; i < 100; i++) {
    		operator.doOperation();
			String newHap = haplotypeModel.getHaplotypeString(0);
			assertNotEquals(haps[0], newHap);
			assertTrue( newHap.contains("C") || newHap.contains("G"));
//			haps[0] = newHap;
			haplotypeModel.reject();
			newHap = haplotypeModel.getHaplotypeString(0);
			assertEquals(haps[0], newHap);

		}

	}
	@Test
	public void testDoOperationMCMC() {
		String[] seqs = new String[]{
				"AAAAAAAAATGTGTTTT....",
				".....CCCCCCCCCCCCCCCCCCCTTTTCCCC....",
				"..........GGGGGGGGGGGGGGCGCGTATAGGGG",
				"...............TTTTTTTTTACACTATA....",
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG"
//				 AAAAACCCCCGCGCCTTCGGTCGTTTTCTATAGGGG
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		String[] haps = new String[]{
				"AAAAACCCCCGGGGGTTTTTACGTACACTATATATA",
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
//				"AAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, hapAlignment);

		
    	// Operators
    	OperatorSchedule schedule = new SimpleOperatorSchedule();

    	MCMCOperator operator = new BasesMultiEmpiricalOperator(haplotypeModel, 3, CoercionMode.COERCION_OFF);
    	operator.setWeight(3.0);
    	schedule.addOperator(operator);
    	
    	//CompoundLikelihood
    	List<Likelihood> likelihoods = new ArrayList<Likelihood>();        

        ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(haplotypeModel);
    	likelihoods.add(srpLikelihood);
    	Likelihood shortReadlikelihood = new CompoundLikelihood(-1, likelihoods);
//    	
//    	
//        likelihoods.clear();
//        likelihoods.add(srpLikelihood);
//        Likelihood posterior = new CompoundLikelihood(0, likelihoods);
//        posterior.setId(CompoundLikelihoodParser.POSTERIOR);
//        
//    	
		double expectedInit = shortReadlikelihood.getLogLikelihood();
		assertEquals(expectedInit, srpLikelihood.getLogLikelihood(), 0);
  	
		int lengthScaler = 1;
    	
		ArrayLogFormatter formatter = new ArrayLogFormatter(false);
    	MCLogger[] loggers = new MCLogger[1];
    	loggers[0] = new MCLogger(formatter, lengthScaler*1, false);
    	loggers[0].add(srpLikelihood );

    	// MCMC
    	MCMC mcmc = new MCMC("mcmc1");
    	MCMCOptions options = new MCMCOptions(lengthScaler*100);
//    	options.setChainLength(10000);
//    	options.setChainLength(lengthScaler*100);
//    	options.setUseCoercion(true); // autoOptimize = true
//    	options.setCoercionDelay(lengthScaler*5);
//    	options.setTemperature(1.0);
//    	options.setFullEvaluationCount(lengthScaler*2);
    	mcmc.setShowOperatorAnalysis(true);
    	mcmc.init(options, srpLikelihood, schedule, loggers);
    	mcmc.run();

    	// time
//    	System.out.println(mcmc.getTimer().toString());

    	// Tracer
    	List<Trace> traces = formatter.getTraces();
//    	ArrayTraceList traceList = new ArrayTraceList("RandomLocalClockTest", traces, 0);
//		Trace trace = traces.get(0);

        for (int i = 0; i < schedule.getOperatorCount(); i++) {
        	MCMCOperator op = schedule.getOperator(i);
        	double acceptanceProb = MCMCOperator.Utils.getAcceptanceProbability(op);
        	assertTrue("0 AcceptanceProb",acceptanceProb>0.01);
        }
		for (Trace<?> trace : traces) {
			if (trace.getName().equals("ShortReadLikelihood")) {

				double startValue = (Double) trace.getValue(0);
				double endValue = (Double) trace
						.getValue(trace.getValuesSize() - 1);
				assertEquals(expectedInit , startValue,0);
				assertTrue(endValue > startValue);
			}
		}

	}



}
