package test.haplotypes.operator;

import static org.junit.Assert.*;


import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;


import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.operator.AlignmentSwapBaseOperator;
import srp.likelihood.ShortReadLikelihood;


import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.tree.Tree.MissingTaxonException;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.TreeLogger;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.ext.SitePatternsExt;
import dr.ext.TreeLikelihoodExt;
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
import dr.inference.operators.Scalable;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.operators.UniformOperator;
import dr.inference.operators.UpDownOperator;
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
    	operator.setWeight(3.0);
    	schedule.addOperator(operator);
    	

    	//CompoundLikelihood

    	
    	List<Likelihood> likelihoods = new ArrayList<Likelihood>();        
//        likelihoods.add(coalescent);
//        Likelihood prior = new CompoundLikelihood(0, likelihoods);
//        prior.setId(CompoundLikelihoodParser.PRIOR);

//        likelihoods.clear();
//        likelihoods.add(treeLikelihood);
//        Likelihood likelihood = new CompoundLikelihood(-1, likelihoods);

//        likelihoods.clear();
        ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(aMap, haplotypeModel);
    	likelihoods.add(srpLikelihood);
    	Likelihood shortReadlikelihood = new CompoundLikelihood(-1, likelihoods);
    	
    	
        likelihoods.clear();
//        likelihoods.add(prior);
//        likelihoods.add(likelihood);
        likelihoods.add(srpLikelihood);
        
        
        Likelihood posterior = new CompoundLikelihood(0, likelihoods);
        posterior.setId(CompoundLikelihoodParser.POSTERIOR);
        
    	
		double expectedInit = shortReadlikelihood.getLogLikelihood();
		assertEquals(expectedInit, srpLikelihood.getLogLikelihood(), 0);

    	ArrayLogFormatter formatter = new ArrayLogFormatter(false);
    	
    	int lengthScaler = 1;
    	MCLogger[] loggers = new MCLogger[1];
    	loggers[0] = new MCLogger(formatter, lengthScaler*1, false);
    	loggers[0].add(shortReadlikelihood );
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



	@Test
	public void testFullMCMCOperation() throws MissingTaxonException, IOException {
		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/";
		
		String trueAlignmentFile = "130220_H4_haplotypes.phyml";
		String truePhylogenyFile = "130220_H4_haplotypes.tree";
		String shortReadFile = "130220_H4_srp.fasta";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		Alignment alignment = dataImporter.importAlignment(trueAlignmentFile);
	
		Tree truePhylogeny = dataImporter.importTree(truePhylogenyFile);
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);
		
		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
		
		
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, alignment);
		{int max = haplotypeModel.getHaplotypeCount() * haplotypeModel.getHaplotypeLength() *1;
		System.out.println("Max swap: "+max);
		for (int i = 0; i < max; i++) {
			haplotypeModel.swapBase();
		}}
		
		
		Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 380.0, 0, 38000.0);
		ConstantPopulationModel startingTree = new ConstantPopulationModel(popSize, Units.Type.YEARS);
		CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel,null, new ArrayList<TaxonList>(), startingTree);
		coalescent.setId("coalescent");

	        // clock model
        Parameter rateParameter =  new Parameter.Default(StrictClockBranchRates.RATE, 1E-5, 0, 100.0);
        StrictClockBranchRates branchRateModel = new StrictClockBranchRates(rateParameter);

	        // Sub model
        Parameter freqs = new Parameter.Default(haplotypeModel.getStateFrequencies());
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100.0);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(hky);
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteModel.setMutationRateParameter(mu);

        //treeLikelihood
        SitePatternsExt patterns = new SitePatternsExt(haplotypeModel, null, 0, -1, 1, true);

        TreeLikelihoodExt treeLikelihood = new TreeLikelihoodExt(patterns, treeModel, siteModel, branchRateModel, null,
                false, false, true, false, false);
        treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);
        
        
    	// Operators
		OperatorSchedule schedule = new SimpleOperatorSchedule();
		schedule = defalutOperators(schedule, kappa, rateParameter, treeModel);
		
    	int index = 0;
    	MCMCOperator operator = new AlignmentSwapBaseOperator(haplotypeModel, index, CoercionMode.COERCION_ON);
    	operator.setWeight(10);
    	schedule.addOperator(operator);
    	
    	//CompoundLikelihood
    	List<Likelihood> likelihoods = new ArrayList<Likelihood>();        
        likelihoods.add(coalescent);
        Likelihood prior = new CompoundLikelihood(0, likelihoods);
        prior.setId(CompoundLikelihoodParser.PRIOR);

        likelihoods.clear();
        likelihoods.add(treeLikelihood);
        Likelihood likelihood = new CompoundLikelihood(-1, likelihoods);

        likelihoods.clear();
        ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(aMap, haplotypeModel);
    	likelihoods.add(srpLikelihood);
    	Likelihood shortReadlikelihood = new CompoundLikelihood(-1, likelihoods);
    	
    	
        likelihoods.clear();
        likelihoods.add(prior);
        likelihoods.add(likelihood);
        likelihoods.add(srpLikelihood);
        
        
        Likelihood posterior = new CompoundLikelihood(0, likelihoods);
        posterior.setId(CompoundLikelihoodParser.POSTERIOR);
        
    	
		double expectedInit = shortReadlikelihood.getLogLikelihood();
		assertEquals(expectedInit, srpLikelihood.getLogLikelihood(), 0);

    	ArrayLogFormatter formatter = new ArrayLogFormatter(false);
    	
    	int lengthScaler = 1000;
    	MCLogger[] loggers = new MCLogger[3];
    	loggers[0] = new MCLogger(formatter, lengthScaler*1, false);
//    	loggers[0] = new MCLogger(new TabDelimitedFormatter(System.out), lengthScaler*1, false);
//    	loggers[0] = new MCLogger("/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/zzzout.log", lengthScaler*1, false, 0);
    	loggers[0].add(prior);
    	loggers[0].add(likelihood);
    	loggers[0].add(shortReadlikelihood );
    	loggers[0].add(posterior);
//    	loggers[0].add(kappa);


    	loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), lengthScaler*1, false);
//    	loggers[1] = new MCLogger("/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/zzzout.log", lengthScaler*1, false, 0);
    	loggers[1].add(prior);
    	loggers[1].add(likelihood);
    	loggers[1].add(shortReadlikelihood );
    	loggers[1].add(posterior);


        File file = new File("testSwap.trees");
        
        final PrintWriter pw = new PrintWriter(new FileOutputStream(file));

        TabDelimitedFormatter treeFormatter = new TabDelimitedFormatter(pw);
        
//        loggers[1] = new TreeLogger(treeModel, new TabDelimitedFormatter(out), lengthScaler, true, true, false);
        loggers[2] = new TreeLogger(treeModel, branchRateModel, null, null, treeFormatter, lengthScaler, true, true, true, null, null/*, Double.NaN*/);

//    	loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), lengthScaler*1, false);
//    	loggers[1].add(likelihood);
//    	loggers[1].add(srpLikelihood);
//    	loggers[1].add(posterior);
//    	loggers[1].add(kappa);

    	// MCMC
    	MCMC mcmc = new MCMC("mcmc1");
    	MCMCOptions options = new MCMCOptions();
//    	options.setChainLength(10000);
    	options.setChainLength(lengthScaler*100);
//    	options.setUseCoercion(true); // autoOptimize = true
//    	options.setCoercionDelay(lengthScaler*5);
//    	options.setTemperature(1.0);
//    	options.setFullEvaluationCount(lengthScaler*2);
    	mcmc.setShowOperatorAnalysis(true);
    	mcmc.init(options, posterior, schedule, loggers);
    	mcmc.run();
    	
//    	out.flush();
//    	out.close();
    	// time
    	System.out.println(mcmc.getTimer().toString());

    	// Tracer
//    	List<Trace> traces = formatter.getTraces();
//    	ArrayTraceList traceList = new ArrayTraceList("test", traces, 0);
//
//    	
//    	
////		Trace trace = traces.get(0);
//		for (Trace trace : traces) {
//			if (trace.getName().equals("ShortReadLikelihood")) {
//
//				double startValue = (Double) trace.getValue(0);
//				double endValue = (Double) trace
//						.getValue(trace.getValuesSize() - 1);
//				assertEquals(expectedInit , startValue,0);
//				assertTrue(endValue > startValue);
////				System.out.println(trace.getName());
////				break;
//			}
//		}
		
//		for (int j = 0; j < trace.getValuesSize(); j++) {
//			System.out.print(trace.getValue(j) +"\t");
//		}
//		System.out.println();
//			System.out.println(Arrays.toString(trace.getRange()));
//			System.out.println(trace.getTraceT9ype());
			
	}

	private static OperatorSchedule defalutOperators(OperatorSchedule schedule, Parameter kappa,
			Parameter rateParameter, TreeModel treeModel) {
		
		
//		MCMCOperator operator = new ScaleOperator(kappa, 0.75);
//		schedule.addOperator(operator);
//		
//		operator = new ScaleOperator(rateParameter, 0.75);
//		operator.setWeight(3.0);
//		schedule.addOperator(operator);
//
//		Parameter allInternalHeights = treeModel.createNodeHeightsParameter(true, true, false);
//		operator = new UpDownOperator(new Scalable[] { new Scalable.Default(
//				rateParameter) }, new Scalable[] { new Scalable.Default(
//				allInternalHeights) }, 0.75, 3.0, CoercionMode.COERCION_ON);
//		schedule.addOperator(operator);
//
//		// operator = new ScaleOperator(popSize, 0.75);
//		// operator.setWeight(3.0);
//		// schedule.addOperator(operator);
//
//		Parameter rootHeight = treeModel.getRootHeightParameter();
//		rootHeight.setId("TREE_HEIGHT");
//		operator = new ScaleOperator(rootHeight, 0.75);
//		operator.setWeight(3.0);
//		schedule.addOperator(operator);
//
//		Parameter internalHeights = treeModel.createNodeHeightsParameter(false,true, false);
//		operator = new UniformOperator(internalHeights, 30.0);
//		schedule.addOperator(operator);
//
//		operator = new SubtreeSlideOperator(treeModel, 15.0, 1.0, true, false,false, false, CoercionMode.COERCION_ON);
//		schedule.addOperator(operator);
//
//		operator = new ExchangeOperator(ExchangeOperator.NARROW, treeModel,15.0);
//		schedule.addOperator(operator);
//
//		operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 3.0);
//		schedule.addOperator(operator);
//
//		operator = new WilsonBalding(treeModel, 3.0);
//		schedule.addOperator(operator);

		return schedule;
	}

}
