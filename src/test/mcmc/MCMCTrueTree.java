package test.mcmc;

import static org.junit.Assert.assertEquals;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.HaplotypeLogger;
import srp.haplotypes.HaplotypeLoggerWithTrueHaplotype;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.operator.HaplotypeRecombinationOperator;
import srp.haplotypes.operator.HaplotypeSwapSectionOperator;
import srp.haplotypes.operator.SingleBaseOperator;
import srp.haplotypes.operator.SwapBasesEmpiricalOperator;
import srp.haplotypes.operator.SwapBasesMultiOperator;
import srp.haplotypes.operator.SingleBaseUniformOperator;
import srp.haplotypes.operator.SwapBasesUniformOperator;
import srp.likelihood.ShortReadLikelihood;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.TreeLogger;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.ext.SitePatternsExt;
import dr.ext.TreeLikelihoodExt;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.OneOnXPrior;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.math.distributions.LogNormalDistribution;
import dr.util.Attribute;

public class MCMCTrueTree {

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
	public void testMCMCFixTree() throws Exception {
		String dataDir = "/home/sw167/workspace/ABI/data/";
		String truePhylogenyFile = "H6_005_true_tree.trees";
		String shortReadFile = "H6_srp.fasta";
		
		String prefix = "FixTree_";
		String logTracerName = prefix+"H6_300_tracer.log";
		String logTreeName = prefix+"H6_300.trees";
		String logHaplotypeName = prefix+"H6_300_haplatype.log";
		String operatorAnalysisFile = prefix+"operatorAnalysisFile";
		
		int logInterval = 1000;
		
		DataImporter dataImporter = new DataImporter(dataDir);
		Tree truePhylogeny = dataImporter.importTree(truePhylogenyFile);
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);
		
		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, 6);

		// coalescent
		Parameter popSize = new Parameter.Default(
				ConstantPopulationModelParser.POPULATION_SIZE, 3000.0, 100, 100000.0);
		ConstantPopulationModel startingTree = new ConstantPopulationModel(popSize, Units.Type.YEARS);

		CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel,null, new ArrayList<TaxonList>(), startingTree);
		coalescent.setId("coalescent");


		
		// clock model
		Parameter rateParameter = new Parameter.Default(StrictClockBranchRates.RATE, 1E-5, 0, 1);
		StrictClockBranchRates branchRateModel = new StrictClockBranchRates(rateParameter);
		// sub model
		Parameter freqs = new Parameter.Default("Frequency", haplotypeModel.getStateFrequencies());
		
		// treeLikelihood
		Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100.0);
		TreeLikelihoodExt treeLikelihood = setupTreeLikelihood(kappa, freqs,
				haplotypeModel, treeModel, branchRateModel);

		// ShortReadLikelihood
		ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(haplotypeModel);

		// CompoundLikelihood
		HashMap<String, Likelihood> compoundlikelihoods = setupCompoundLikelihood(
				popSize, kappa, coalescent, treeLikelihood, srpLikelihood);
		
		Likelihood prior = compoundlikelihoods.get(CompoundLikelihoodParser.PRIOR);
		Likelihood likelihood = compoundlikelihoods.get(CompoundLikelihoodParser.LIKELIHOOD);
		Likelihood shortReadLikelihood = compoundlikelihoods.get(ShortReadLikelihood.SHORT_READ_LIKELIHOOD);
		Likelihood posterior = compoundlikelihoods.get(CompoundLikelihoodParser.POSTERIOR);
		
		// Operators
		OperatorSchedule schedule = new SimpleOperatorSchedule();
		schedule = defalutOperators(schedule, haplotypeModel, popSize, kappa);
		Parameter rootHeight = treeModel.getRootHeightParameter();
		
		
		double expectedInit = shortReadLikelihood.getLogLikelihood();
		assertEquals(expectedInit, srpLikelihood.getLogLikelihood(), 0);

		
		
//		ArrayLogFormatter formatter = new ArrayLogFormatter(false);
		MCLogger[] loggers = new MCLogger[4];
//		loggers[0] = new MCLogger(formatter, logInterval, false);
		
		loggers[0] = new MCLogger(logTracerName, logInterval, false, 0);
		addToLogger(loggers[0], posterior, prior, likelihood, shortReadLikelihood,
//				rootHeight, rateParameter,
				popSize, kappa, coalescent,
				freqs
				);

		loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), logInterval, true, logInterval*2);
		addToLogger(loggers[1], posterior, prior, likelihood, shortReadLikelihood,
				popSize, kappa,				coalescent);
		
		TabDelimitedFormatter treeFormatter = new TabDelimitedFormatter(
				new PrintWriter(new FileOutputStream(new File(logTreeName))));

		loggers[2] = new TreeLogger(treeModel, branchRateModel, null, null,
				treeFormatter, logInterval, true, true, true, null, null);

		Alignment trueAlignment = dataImporter.importAlignment("H6_005.fasta");
		loggers[3] = new HaplotypeLoggerWithTrueHaplotype(haplotypeModel, trueAlignment, logHaplotypeName, logInterval*10);

		// MCMC
		
		MCMCOptions options = setMCMCOptions(logInterval);
		
		MCMC mcmc = new MCMC("mcmc1");
		mcmc.setShowOperatorAnalysis(true);
		mcmc.setOperatorAnalysisFile(new File(operatorAnalysisFile));
		
		mcmc.init(options, posterior, schedule, loggers);
		mcmc.run();

		System.out.println(mcmc.getTimer().toString());
		// Tracer
		// List<Trace> traces = formatter.getTraces();
		// ArrayTraceList traceList = new ArrayTraceList("test", traces, 0);
		//
		//
		//
		// // Trace trace = traces.get(0);
		// for (Trace trace : traces) {
		// if (trace.getName().equals("ShortReadLikelihood")) {
		//
		// double startValue = (Double) trace.getValue(0);
		// double endValue = (Double) trace
		// .getValue(trace.getValuesSize() - 1);
		// assertEquals(expectedInit , startValue,0);
		// assertTrue(endValue > startValue);
		// // System.out.println(trace.getName());
		// // break;
		// }
		// }

		// for (int j = 0; j < trace.getValuesSize(); j++) {
		// System.out.print(trace.getValue(j) +"\t");
		// }
		// System.out.println();
		// System.out.println(Arrays.toString(trace.getRange()));
		// System.out.println(trace.getTraceT9ype());

	}

	private TreeLikelihoodExt setupTreeLikelihood(Parameter kappa, Parameter freqs,
			HaplotypeModel haplotypeModel, TreeModel treeModel, StrictClockBranchRates branchRateModel) {

		// Sub model
		
		FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
		HKY hky = new HKY(kappa, f);

		// siteModel
		GammaSiteModel siteModel = new GammaSiteModel(hky);
		Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1, 0, Double.POSITIVE_INFINITY);
		siteModel.setMutationRateParameter(mu);

		// treeLikelihood
//		SitePatternsExt patterns = new SitePatternsExt(haplotypeModel, null, 0, -1, 1, true);

		TreeLikelihoodExt treeLikelihood = new TreeLikelihoodExt(haplotypeModel, treeModel, siteModel,
				 branchRateModel, null, false, false, true, false, false);
//		TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel,
//				siteModel, branchRateModel, null, false, false, true, false, false);
		treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);
		
		return treeLikelihood;
	}

	private HashMap<String, Likelihood> setupCompoundLikelihood(
			Parameter popSize, Parameter kappa, Likelihood coalescent,
			Likelihood treeLikelihood, Likelihood srpLikelihood) {

		//		 CompoundLikelihood
		HashMap<String, Likelihood> compoundLikelihoods = new HashMap<String, Likelihood>(4);
		
		List<Likelihood> likelihoods = new ArrayList<Likelihood>();

		// Prior
		OneOnXPrior oneOnX = new OneOnXPrior();
		oneOnX.addData(popSize);

		DistributionLikelihood logNormalLikelihood = new DistributionLikelihood(
				new LogNormalDistribution(1.0, 1.25), 0); // meanInRealSpace="false"
		logNormalLikelihood.addData(kappa);

		likelihoods.add(oneOnX);
		likelihoods.add(logNormalLikelihood);
		likelihoods.add(coalescent);
		Likelihood prior = new CompoundLikelihood(0, likelihoods);
		prior.setId(CompoundLikelihoodParser.PRIOR);
		compoundLikelihoods.put(CompoundLikelihoodParser.PRIOR, prior);
		
		// Likelihood
		likelihoods.clear();
		likelihoods.add(treeLikelihood);
		Likelihood likelihood = new CompoundLikelihood(-1, likelihoods);
		likelihood.setId(CompoundLikelihoodParser.LIKELIHOOD);
		compoundLikelihoods.put(CompoundLikelihoodParser.LIKELIHOOD, likelihood);
		
		// ShortReadLikelihood
		likelihoods.clear();
		
		likelihoods.add(srpLikelihood);
		Likelihood shortReadlikelihood = new CompoundLikelihood(-1, likelihoods);
		shortReadlikelihood.setId(ShortReadLikelihood.SHORT_READ_LIKELIHOOD);
		compoundLikelihoods.put(ShortReadLikelihood.SHORT_READ_LIKELIHOOD, shortReadlikelihood);
		
		// Posterior
		likelihoods.clear();
		likelihoods.add(prior);
		likelihoods.add(likelihood);
		likelihoods.add(srpLikelihood);
		Likelihood posterior = new CompoundLikelihood(0, likelihoods);
		posterior.setId(CompoundLikelihoodParser.POSTERIOR);
		compoundLikelihoods.put(CompoundLikelihoodParser.POSTERIOR, posterior);

		return compoundLikelihoods;
	}

	private static MCMCOptions setMCMCOptions(int logInterval) {
		MCMCOptions options = new MCMCOptions();
		options.setChainLength(logInterval * 10);
		options.setUseCoercion(true); // autoOptimize = true
		options.setCoercionDelay(logInterval * 0);
		options.setTemperature(1.0);
		options.setFullEvaluationCount(logInterval*0);

		return options;
	}

	private static OperatorSchedule defalutOperators(OperatorSchedule schedule,
			HaplotypeModel haplotypeModel, Parameter popSize, Parameter kappa) {

		MCMCOperator operator;
		
		operator = new SingleBaseOperator(haplotypeModel,0);
		operator.setWeight(3.0);
		schedule.addOperator(operator);

//		operator = new SingleBaseUniformOperator(haplotypeModel,0);
//		operator.setWeight(3.0);
//		schedule.addOperator(operator);

		operator = new SwapBasesMultiOperator(haplotypeModel, 2, CoercionMode.COERCION_OFF);
		operator.setWeight(3.0); 
		schedule.addOperator(operator);


//		operator = new SwapBasesUniformOperator(haplotypeModel, 1, CoercionMode.COERCION_OFF);
//		operator.setWeight(3.0); 
//		schedule.addOperator(operator);
//
//		operator = new SwapBasesEmpiricalOperator(haplotypeModel, 6, CoercionMode.COERCION_ON);
//		operator.setWeight(6.0); 
//		schedule.addOperator(operator);
//		
		operator = new HaplotypeRecombinationOperator(haplotypeModel, 12);
		operator.setWeight(6.0); 
		schedule.addOperator(operator);
//		
//		operator = new HaplotypeSwapSectionOperator(haplotypeModel, 24, CoercionMode.COERCION_ON);
//		operator.setWeight(6.0); 
//		schedule.addOperator(operator);
//		
		
//		operator = new ScaleOperator(kappa, 0.75);
//		operator.setWeight(0.1);
//		schedule.addOperator(operator);

//		operator = new ScaleOperator(popSize, 0.75);
//		operator.setWeight(0.1);
//		schedule.addOperator(operator);


		return schedule;
	}

	public static MCLogger addToLogger(MCLogger mcLogger, Loggable... loggableParameter) {
		for (Loggable loggable : loggableParameter) {
			mcLogger.add(loggable);
		}
		return mcLogger;
		
	}

	}
