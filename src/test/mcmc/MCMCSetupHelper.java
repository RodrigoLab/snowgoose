package test.mcmc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.operator.HaplotypeRecombinationOperator;
import srp.haplotypes.operator.HaplotypeSwapSectionOperator;
import srp.haplotypes.operator.SingleBaseEmpiricalOperator;
import srp.haplotypes.operator.SingleBaseFrequencyOperator;
import srp.haplotypes.operator.SingleBaseOperator;
import srp.haplotypes.operator.SingleBaseUniformOperator;
import srp.haplotypes.operator.SwapBasesEmpiricalOperator;
import srp.haplotypes.operator.SwapBasesMultiOperator;
import srp.haplotypes.operator.SwapBasesUniformOperator;
import srp.likelihood.ShortReadLikelihood;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evolution.util.Units.Type;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.ext.TreeLikelihoodExt;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.MCLogger;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.OneOnXPrior;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.DeltaExchangeOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.Scalable;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.UniformOperator;
import dr.inference.operators.UpDownOperator;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.math.distributions.LogNormalDistribution;

public class MCMCSetupHelper {

	public static TreeLikelihoodExt setupTreeLikelihood(Parameter kappa,
			Parameter freqs, HaplotypeModel haplotypeModel,
			TreeModel treeModel, StrictClockBranchRates branchRateModel) {

		// Sub model
		FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
		HKY hky = new HKY(kappa, f);

		// siteModel
		GammaSiteModel siteModel = new GammaSiteModel(hky);
		Parameter mu = new Parameter.Default(
				GammaSiteModelParser.MUTATION_RATE, 1, 0, Double.POSITIVE_INFINITY);
		siteModel.setMutationRateParameter(mu);

		// treeLikelihood

		TreeLikelihoodExt treeLikelihood = new TreeLikelihoodExt(
				haplotypeModel, treeModel, siteModel, branchRateModel, null,
				false, false, true, false, false);
		treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);
		// SitePatternsExt patterns = new SitePatternsExt(haplotypeModel, null, 0, -1, 1, true);
		// TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null, false, false, true, false, false);
		

		return treeLikelihood;
	}

	public static HashMap<String, Likelihood> setupCompoundLikelihood(
			Parameter popSize, Parameter kappa, Likelihood coalescent,
			Likelihood treeLikelihood, Likelihood srpLikelihood) {

		// CompoundLikelihood
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
		compoundLikelihoods.put(ShortReadLikelihood.SHORT_READ_LIKELIHOOD,shortReadlikelihood);

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
	private static final double opTiny = 0.1;
	private static final double opSmall = 3;
	private static final double opMed = 10;
	private static final double opLarge = 30;
//	private static final double opHuge = 300;
	
	public static ArrayList<MCMCOperator> defalutOperators(HaplotypeModel haplotypeModel,
			Parameter... parameters) {

		MCMCOperator operator;
		ArrayList<MCMCOperator> OperatorList = new ArrayList<MCMCOperator>();

//		operator = new SingleBaseOperator(haplotypeModel, 0);
//		operator.setWeight(opSmall);
//		OperatorList.add(operator);
//
//		operator = new SingleBaseUniformOperator(haplotypeModel, 0);
//		operator.setWeight(opSmall);
//		OperatorList.add(operator);

		operator = new SingleBaseEmpiricalOperator(haplotypeModel, 0);
		operator.setWeight(opLarge);
		OperatorList.add(operator);

//		operator = new SwapBasesMultiOperator(haplotypeModel, 12, CoercionMode.COERCION_ON);
//		operator.setWeight(opLarge);
//		OperatorList.add(operator);
//
//		operator = new SwapBasesUniformOperator(haplotypeModel, 12, CoercionMode.COERCION_ON);
//		operator.setWeight(opLarge);
//		OperatorList.add(operator);
//
		operator = new SwapBasesEmpiricalOperator(haplotypeModel, 3, CoercionMode.COERCION_ON);
		operator.setWeight(opLarge);
		OperatorList.add(operator);

//		operator = new HaplotypeRecombinationOperator(haplotypeModel, 12);
//		operator.setWeight(3.0);
//		OperatorList.add(operator);

//		operator = new HaplotypeSwapSectionOperator(haplotypeModel, 12, CoercionMode.COERCION_ON);
//		operator.setWeight(opSmall);
//		OperatorList.add(operator);
		
		for (Parameter parameter : parameters) {
			String parameterName = parameter.getParameterName();
			
			if("kappa".equals(parameterName)){
				operator = new ScaleOperator(parameter, 0.75);
				operator.setWeight(opTiny);
				OperatorList.add(operator);
			}
			else if("frequency".equals(parameterName)){
				operator = new DeltaExchangeOperator(parameter, new int[] { 1,
						1, 1, 1 }, 0.01, 0.1, false, CoercionMode.COERCION_ON);
				operator.setWeight(opTiny);
				OperatorList.add(operator);
				

				operator = new SingleBaseFrequencyOperator(haplotypeModel, parameter);
				operator.setWeight(opLarge);
				OperatorList.add(operator);
			}
			else if("populationSize".equals(parameterName)){
				operator = new ScaleOperator(parameter, 0.75);
				operator.setWeight(opTiny);
				OperatorList.add(operator);
			}
			
		}
		
		return OperatorList;
	}

	
	public static List<MCMCOperator> defalutTreeOperators(TreeModel treeModel) {

		MCMCOperator operator;
		List<MCMCOperator> OperatorList = new ArrayList<MCMCOperator>();
		
		Parameter allInternalHeights = treeModel.createNodeHeightsParameter(true, true, false);
		operator = new UpDownOperator(
				null,// new Scalable[] { new Scalable.Default(rateParameter) },
				new Scalable[] { new Scalable.Default(allInternalHeights) },
				0.75, 3.0, CoercionMode.COERCION_ON);
		operator.setWeight(opSmall);
		OperatorList.add(operator);

		Parameter rootHeight = treeModel.getRootHeightParameter();
		rootHeight.setId("TREE_HEIGHT");
		operator = new ScaleOperator(rootHeight, 0.75);
		operator.setWeight(opSmall);
		OperatorList.add(operator);

		Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);
		operator = new UniformOperator(internalHeights, 30.0);
		operator.setWeight(opMed);
		OperatorList.add(operator);

		operator = new SubtreeSlideOperator(treeModel, 15.0, 1.0, true, false,
				false, false, CoercionMode.COERCION_ON);
		operator.setWeight(opMed);
		OperatorList.add(operator);

		operator = new ExchangeOperator(ExchangeOperator.NARROW, treeModel, 15.0);
		operator.setWeight(opMed);
		OperatorList.add(operator);

		operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 3.0);
		operator.setWeight(opSmall);
		OperatorList.add(operator);

		operator = new WilsonBalding(treeModel, 3.0);
		operator.setWeight(opSmall);
		OperatorList.add(operator);

		return OperatorList;
	}

	public static MCLogger addToLogger(MCLogger mcLogger,
			Loggable... loggableParameter) {
		for (Loggable loggable : loggableParameter) {
			mcLogger.add(loggable);
		}
		return mcLogger;

	}

	public static TreeModel setupRandomTreeModel(ConstantPopulationModel startingTree, HaplotypeModel haplotypeModel, Type years) {
		
		ConstantPopulation constant = (ConstantPopulation) startingTree.getDemographicFunction();
		CoalescentSimulator simulator = new CoalescentSimulator();
		Tree tree = simulator.simulateTree(haplotypeModel, constant);
		TreeModel treeModel = new TreeModel(tree);// treeModel

		return treeModel;
	}

	public static MCMCOptions setMCMCOptions(int logInterval, int totalSamples) {
		MCMCOptions options = new MCMCOptions();
		options.setChainLength(logInterval * totalSamples);
		options.setUseCoercion(false); // autoOptimize = true
		options.setCoercionDelay(logInterval * 2);
		options.setTemperature(1.0);
		options.setFullEvaluationCount(logInterval*0);
	
		return options;
	}


	public static ArrayList<MCMCOperator> testOperators(HaplotypeModel haplotypeModel,
			Parameter... parameters) {

		MCMCOperator operator;
		ArrayList<MCMCOperator> OperatorList = new ArrayList<MCMCOperator>();

//		operator = new SingleBaseOperator(haplotypeModel, 0);
//		operator.setWeight(opSmall);
//		OperatorList.add(operator);
//
//		operator = new SingleBaseUniformOperator(haplotypeModel, 0);
//		operator.setWeight(opSmall);
//		OperatorList.add(operator);

		operator = new SingleBaseEmpiricalOperator(haplotypeModel, 0);
		operator.setWeight(opMed);
		OperatorList.add(operator);

//		operator = new SwapBasesMultiOperator(haplotypeModel, 12, CoercionMode.COERCION_ON);
//		operator.setWeight(opLarge);
//		OperatorList.add(operator);
//
//		operator = new SwapBasesUniformOperator(haplotypeModel, 12, CoercionMode.COERCION_ON);
//		operator.setWeight(opLarge);
//		OperatorList.add(operator);
//
//		operator = new SwapBasesEmpiricalOperator(haplotypeModel, 3, CoercionMode.COERCION_ON);
//		operator.setWeight(opLarge);
//		OperatorList.add(operator);

//		operator = new HaplotypeRecombinationOperator(haplotypeModel, 12);
//		operator.setWeight(3.0);
//		OperatorList.add(operator);

//		operator = new HaplotypeSwapSectionOperator(haplotypeModel, 12, CoercionMode.COERCION_ON);
//		operator.setWeight(opSmall);
//		OperatorList.add(operator);
		
		for (Parameter parameter : parameters) {
			String parameterName = parameter.getParameterName();
			
			if("kappa".equals(parameterName)){
				operator = new ScaleOperator(parameter, 0.75);
				operator.setWeight(opTiny);
//				OperatorList.add(operator);
			}
			else if("frequency".equals(parameterName)){
				operator = new DeltaExchangeOperator(parameter, new int[] { 1,
						1, 1, 1 }, 0.01, 0.1, false, CoercionMode.COERCION_ON);
				operator.setWeight(opTiny);
//				operator.setWeight(opLarge);
				OperatorList.add(operator);
				

				operator = new SingleBaseFrequencyOperator(haplotypeModel, parameter);
				operator.setWeight(opMed);
				OperatorList.add(operator);
			}
			else if("populationSize".equals(parameterName)){
				operator = new ScaleOperator(parameter, 0.75);
				operator.setWeight(opTiny);
				OperatorList.add(operator);
			}
			
		}
		
		return OperatorList;
	}

	
}
