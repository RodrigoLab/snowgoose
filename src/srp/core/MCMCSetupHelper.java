package srp.core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import srp.likelihood.ShortReadLikelihood;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units.Type;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.tree.TreeModel;
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
		
			
			MCMCOptions options = new MCMCOptions(logInterval * totalSamples, 10,
					1, MarkovChain.EVALUATION_TEST_THRESHOLD, true, 0, 1.0);
			//		MCMCOptions(long chainLength, 
	//				long fullEvaluationCount, 
	//				int minOperatorCountForFullEvaluation, 
	//				double evaluationTestThreshold, 
	//				boolean coercion, 
	//				long coercionDelay, 
	//				double temperature) 
			return options;
		    
		}

	public MCMCSetupHelper() {
		super();
	}

}