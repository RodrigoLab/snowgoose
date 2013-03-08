package test.srp.haplotypes.operator;

import static org.junit.Assert.*;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.operator.AlignmentSwapBaseOperator;
import srp.likelihood.ShortReadLikelihood;


import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evomodelxml.substmodel.HKYParser;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.trace.ArrayTraceList;
import dr.inference.trace.Trace;
import dr.inferencexml.model.CompoundLikelihoodParser;

public class AlignmentSwapBaseOperatorTest {

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
//		fail("Not yet implemented");
	}

	@Test
	public void testDoOperation() {
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
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, hapAlignment);

		
    	// Operators
    	OperatorSchedule schedule = new SimpleOperatorSchedule();

    	Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100.0);
    	MCMCOperator operator = new ScaleOperator(kappa, 0.75);
    	schedule.addOperator(operator);

    	int index = 0;
    	
    	operator = new AlignmentSwapBaseOperator(haplotypeModel, index, CoercionMode.COERCION_ON);
    	operator.setWeight(10.0);
    	schedule.addOperator(operator);
    	

    	//CompoundLikelihood
    	List<Likelihood> likelihoodsList = new ArrayList<Likelihood>();       
		ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(aMap, haplotypeModel);
//    	likelihoods.add(treeLikelihood);
    	likelihoodsList.add(srpLikelihood);

    	Likelihood likelihood = new CompoundLikelihood(0, likelihoodsList);
    	likelihood.setId("ShortReadLikelihood");
    	
    	Likelihood posterior = new CompoundLikelihood(0, likelihoodsList);
    	posterior.setId(CompoundLikelihoodParser.POSTERIOR);

		double expectedInit = likelihood.getLogLikelihood();
		assertEquals(expectedInit, srpLikelihood.getLogLikelihood(), 0);
    	//    	likelihoods.clear();
//    	likelihoods.add(prior);
//    	likelihoods.add(likelihood);
//    	likelihoods.add(shortReadLikelihood);
//    	Likelihood posterior = new CompoundLikelihood(0, likelihoods);
//    	posterior.setId(CompoundLikelihoodParser.POSTERIOR);

    	ArrayLogFormatter formatter = new ArrayLogFormatter(false);
    	
    	int lengthScaler = 1;
    	MCLogger[] loggers = new MCLogger[1];
    	loggers[0] = new MCLogger(formatter, lengthScaler*1, false);
    	loggers[0].add(likelihood);
    	loggers[0].add(srpLikelihood);
    	loggers[0].add(posterior);
    	loggers[0].add(kappa);

//    	loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), lengthScaler*1, false);
//    	loggers[1].add(likelihood);
//    	loggers[1].add(srpLikelihood);
//    	loggers[1].add(posterior);
//    	loggers[1].add(kappa);

    	// MCMC
    	MCMC mcmc = new MCMC("mcmc1");
    	MCMCOptions options = new MCMCOptions();
//    	options.setChainLength(10000);
    	options.setChainLength(lengthScaler*10);
//    	options.setUseCoercion(true); // autoOptimize = true
//    	options.setCoercionDelay(lengthScaler*5);
//    	options.setTemperature(1.0);
//    	options.setFullEvaluationCount(lengthScaler*2);
//    	mcmc.setShowOperatorAnalysis(true);
    	mcmc.init(options, posterior, schedule, loggers);
    	mcmc.run();

    	// time
//    	System.out.println(mcmc.getTimer().toString());

    	// Tracer
    	List<Trace> traces = formatter.getTraces();
    	ArrayTraceList traceList = new ArrayTraceList("RandomLocalClockTest", traces, 0);

    	
    	
//		Trace trace = traces.get(0);
		for (Trace trace : traces) {
			if (trace.getName().equals("ShortReadLikelihood")) {

				double startValue = (Double) trace.getValue(0);
				double endValue = (Double) trace
						.getValue(trace.getValuesSize() - 1);
				assertEquals(expectedInit , startValue,0);
				assertTrue(endValue > startValue);
//				System.out.println(trace.getName());
//				break;
			}
		}
		
//		for (int j = 0; j < trace.getValuesSize(); j++) {
//			System.out.print(trace.getValue(j) +"\t");
//		}
//		System.out.println();
//			System.out.println(Arrays.toString(trace.getRange()));
//			System.out.println(trace.getTraceT9ype());
			
	}


	

}
