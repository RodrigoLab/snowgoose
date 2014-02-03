package srp.core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import srp.spectrum.likelihood.ShortReadsSpectrumLikelihood;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units.Type;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.MCLogger;
import dr.inference.markovchain.MarkovChain;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.OneOnXPrior;
import dr.inference.model.Parameter;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.math.distributions.LogNormalDistribution;

public class MCMCSetupHelper {

	public static TreeModel setupRandomTreeModel(ConstantPopulationModel popModel, TaxonList taxonList,
			Type years) {
		
		ConstantPopulation constant = (ConstantPopulation) popModel.getDemographicFunction();
		CoalescentSimulator simulator = new CoalescentSimulator();
		Tree tree = simulator.simulateTree(taxonList, constant);
		TreeModel treeModel = new TreeModel(tree);// treeModel
	
		return treeModel;
	}

	public static HashMap<String, Likelihood> setupCompoundLikelihood(Parameter popSize, Parameter kappa,
			Likelihood coalescent, Likelihood treeLikelihood, Likelihood srpLikelihood) {
	
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
		shortReadlikelihood.setId(ShortReadsSpectrumLikelihood.SHORT_READ_LIKELIHOOD);
		compoundLikelihoods.put(ShortReadsSpectrumLikelihood.SHORT_READ_LIKELIHOOD, shortReadlikelihood);
	
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

	public static MCLogger addToLogger(MCLogger mcLogger, Loggable... loggableParameter) {
		for (Loggable loggable : loggableParameter) {
			mcLogger.add(loggable);
		}
		return mcLogger;
	
	}

	public static MCMCOptions setMCMCOptions(int logInterval, int totalSamples) {
	//		MCMCOptions options = new MCMCOptions();
	//		options.setChainLength(logInterval * totalSamples);
	//		options.setUseCoercion(false); // autoOptimize = true
	//		options.setCoercionDelay((int) (logInterval * 0.01));
	//		options.setTemperature(1.0);
	//		options.setFullEvaluationCount((int) (logInterval*0.01));
		
		int coercionDelay = logInterval * totalSamples /100;
		MCMCOptions options = new MCMCOptions(logInterval * totalSamples, 100,
				0, MarkovChain.EVALUATION_TEST_THRESHOLD, false, coercionDelay, 1.0);
			//		MCMCOptions(long chainLength, 
	//				long fullEvaluationCount, //2000
	//				int minOperatorCountForFullEvaluation, //1 
	//				double evaluationTestThreshold, 
	//				boolean coercion, 
	//				long coercionDelay, //chainLength/100 
	//				double temperature) 
		return options;
		    
	}
	
	public static SiteModel setupSiteModel(){
		
		// clock model
		Parameter rateParameter = new Parameter.Default(StrictClockBranchRates.RATE, 1e-5, 0, 1);
		StrictClockBranchRates branchRateModel = new StrictClockBranchRates(rateParameter);
	
		Parameter freqs = new Parameter.Default("frequency", new double[]{0.25, 0.25, 0.25, 0.25});
		Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100.0);

		// Sub model
		FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
		HKY hky = new HKY(kappa, f);

		// siteModel
		GammaSiteModel siteModel = new GammaSiteModel(hky);
		Parameter mu = new Parameter.Default(
				GammaSiteModelParser.MUTATION_RATE, 1, 0, Double.POSITIVE_INFINITY);
		siteModel.setMutationRateParameter(mu);
		return siteModel;
	}
	

}