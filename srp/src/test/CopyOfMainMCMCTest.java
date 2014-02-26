package test;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.junit.Test;

import srp.core.DataImporter;
import srp.core.MCMCSetupHelper;
import srp.core.MCMCSetupHelperSpectrum;
import srp.haplotypes.AlignmentMapping;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumAlignmentUtils;
import srp.spectrum.likelihood.ShortReadsSpectrumLikelihood;
import srp.spectrum.treelikelihood.SpectrumTreeLikelihood;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.OneOnXPrior;
import dr.inference.model.Parameter;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.math.distributions.LogNormalDistribution;

public class CopyOfMainMCMCTest {

	@Test
	public void test() throws Exception{


		String dataDir = "/home/sw167/workspaceSrp/ABI/unittest/testData/";
		int runIndex = 1;
		int totalSamples = 10;
		int logInterval = 1;
		int noOfTrueHaplotype = 7;
		int noOfRecoveredHaplotype=7;
		String hapRunIndex = "H"+noOfTrueHaplotype+"_"+runIndex;
		String shortReadFile = hapRunIndex +"_Srp.fasta";
		String trueHaplotypeFile = hapRunIndex +"_Srp_fullHaplotype.fasta";
		
		String prefix = dataDir+"FullTree_"+hapRunIndex;
		String logTracerName = prefix+".log";
		String logTreeName = prefix+".trees";
		String operatorAnalysisFile = prefix+"_operatorAnalysisFile.txt";
		
		DataImporter dataImporter = new DataImporter(dataDir);

		Alignment trueAlignment = dataImporter.importAlignment(trueHaplotypeFile);
		
		Alignment shortReads = dataImporter.importShortReads(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
		int spectrumLength = aMap.getLength();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, noOfRecoveredHaplotype);
		
		// coalescent
		Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000.0, 100, 100000.0);

		// Random treeModel
		ConstantPopulationModel popModel = new ConstantPopulationModel(popSize, Units.Type.YEARS);
		TreeModel treeModel = MCMCSetupHelper.setupRandomTreeModel(popModel, spectrumModel, Units.Type.YEARS);

		// Coalescent likelihood
		CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel,null, new ArrayList<TaxonList>(), popModel);
		coalescent.setId("coalescent");

		// Simulate haplotypes, treeLikelihood

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

//		// Simulate halotypes
//		if(errorRate>0){
//			haplotypeModel.simulateSequence(errorRate, siteModel, hky, treeModel);
//		}
		
		// treeLikelihood
		SpectrumTreeLikelihood treeLikelihood = new SpectrumTreeLikelihood(
				spectrumModel, treeModel, siteModel, branchRateModel, 
				false, false, true, false, false);
		treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);
		
		// ShortReadLikelihood
		ShortReadsSpectrumLikelihood srpLikelihood = new ShortReadsSpectrumLikelihood(spectrumModel, aMap);

		// CompoundLikelihood
		HashMap<String, Likelihood> compoundlikelihoods = MCMCSetupHelper.setupCompoundLikelihood(
				popSize, kappa, coalescent, treeLikelihood, srpLikelihood);
		Likelihood prior = compoundlikelihoods.get(CompoundLikelihoodParser.PRIOR);
		Likelihood likelihood = compoundlikelihoods.get(CompoundLikelihoodParser.LIKELIHOOD);
		Likelihood shortReadLikelihood = compoundlikelihoods.get(ShortReadsSpectrumLikelihood.SHORT_READ_LIKELIHOOD);
		Likelihood posterior = compoundlikelihoods.get(CompoundLikelihoodParser.POSTERIOR);

		// Operators
		OperatorSchedule schedule = new SimpleOperatorSchedule();
		MCMCSetupHelperSpectrum.defalutSpectrumOperators(schedule, spectrumModel);//, freqs, kappa);
//		MCMCSetupHelperSpectrum.defalutTreeOperators(schedule, treeModel);

		// MCLogger
		MCLogger[] loggers = new MCLogger[1];
		// log tracer
		loggers[0] = new MCLogger(new TabDelimitedFormatter(System.out), logInterval, true, logInterval*2);
		MCMCSetupHelper.addToLogger(loggers[0],
				posterior, prior, likelihood, shortReadLikelihood,
				kappa
				);
		
		// MCMC
		MCMCOptions options = MCMCSetupHelper.setMCMCOptions(logInterval, totalSamples);
		
		MCMC mcmc = new MCMC("mcmc1");
		mcmc.setShowOperatorAnalysis(true);
		mcmc.setOperatorAnalysisFile(new File(operatorAnalysisFile));
		
		System.out.println(likelihood.getLogLikelihood());
		mcmc.init(options, likelihood, schedule, loggers);
		mcmc.run();
		
		System.out.println(mcmc.getTimer().toString());
		
		System.out.println(likelihood.getLogLikelihood());
		
		///////////////////////////////////////////////////////////
	
		
		for (int i = 0; i < 10; i++) {
			

			TreeModel newTreeModel  = new TreeModel(treeModel);
			
			SpectrumAlignmentModel newSpectrumModel = SpectrumAlignmentModel.duplicateSpectrumAlignmentModel(spectrumModel);
			ShortReadsSpectrumLikelihood newSrpLikelihood = new ShortReadsSpectrumLikelihood(newSpectrumModel, aMap);

			// clock model
			Parameter newRateParameter = new Parameter.Default(StrictClockBranchRates.RATE, 1e-5, 0, 1);
			StrictClockBranchRates newBranchRateModel = new StrictClockBranchRates(newRateParameter);
		
			Parameter newFreqs = new Parameter.Default("frequency", freqs.getParameterValues());
			Parameter newKappa = new Parameter.Default(HKYParser.KAPPA, kappa.getParameterValue(0), 0, 100.0);

			// Sub model
			FrequencyModel newFM = new FrequencyModel(Nucleotides.INSTANCE, newFreqs);
			HKY newHky = new HKY(newKappa, newFM);
			// siteModel
			GammaSiteModel newSiteModel = new GammaSiteModel(newHky);
			Parameter newMu = new Parameter.Default(
					GammaSiteModelParser.MUTATION_RATE, 1, 0, Double.POSITIVE_INFINITY);
			newSiteModel.setMutationRateParameter(newMu);
			// treeLikelihood

			
			SpectrumTreeLikelihood atreeLikelihood = new SpectrumTreeLikelihood(
					newSpectrumModel, newTreeModel, newSiteModel, newBranchRateModel, 
					false, false, true, false, false);
			System.out.println(atreeLikelihood.getLogLikelihood());
	}
		/////////////////////////////////////////////////////////
		TreeModel newTreeModel  = new TreeModel(treeModel);
		
		SpectrumAlignmentModel newSpectrumModel = SpectrumAlignmentModel.duplicateSpectrumAlignmentModel(spectrumModel);
		ShortReadsSpectrumLikelihood newSrpLikelihood = new ShortReadsSpectrumLikelihood(newSpectrumModel, aMap);

		// clock model
		Parameter newRateParameter = new Parameter.Default(StrictClockBranchRates.RATE, 1e-5, 0, 1);
		StrictClockBranchRates newBranchRateModel = new StrictClockBranchRates(newRateParameter);
	
		Parameter newFreqs = new Parameter.Default("frequency", freqs.getParameterValues());
		Parameter newKappa = new Parameter.Default(HKYParser.KAPPA, kappa.getParameterValue(0), 0, 100.0);

		// Sub model
		FrequencyModel newFM = new FrequencyModel(Nucleotides.INSTANCE, newFreqs);
		HKY newHky = new HKY(newKappa, newFM);
		System.out.println(Arrays.toString(newFM.getFrequencies()));
		// siteModel
		GammaSiteModel newSiteModel = new GammaSiteModel(newHky);
		Parameter newMu = new Parameter.Default(
				GammaSiteModelParser.MUTATION_RATE, 1, 0, Double.POSITIVE_INFINITY);
		newSiteModel.setMutationRateParameter(newMu);
		// treeLikelihood

		
		System.out.println("================");
		treeLikelihood = new SpectrumTreeLikelihood(
				newSpectrumModel, newTreeModel, newSiteModel, newBranchRateModel, 
				false, false, true, false, false);
		System.out.println(treeLikelihood.getLogLikelihood());

		treeLikelihood = new SpectrumTreeLikelihood(
				spectrumModel, treeModel, siteModel, branchRateModel, 
				false, false, true, false, false);
		System.out.println(treeLikelihood.getLogLikelihood());
		
		treeLikelihood = new SpectrumTreeLikelihood(
				spectrumModel, treeModel, siteModel, branchRateModel, 
				false, false, true, false, false);
		System.out.println(treeLikelihood.getLogLikelihood());
		
		treeLikelihood = new SpectrumTreeLikelihood(
				spectrumModel, treeModel, siteModel, branchRateModel, 
				false, false, true, false, false);
		System.out.println(treeLikelihood.getLogLikelihood());
		
		for (int i = 0; i < treeLikelihood.getModelCount(); i++) {
			Model model = treeLikelihood.getModel(i);
			System.out.println(model.getModelName());
			if(model.getModelName().equals("siteModel")){
				newSiteModel = (GammaSiteModel) model;
				System.out.println(Arrays.toString(newSiteModel.getFrequencyModel()
						.getFrequencies()));
			}
			else if(model.getModelName().equals("treeModel")){
				System.out.println(newTreeModel.toString());
				newTreeModel = (TreeModel) model;
				System.out.println(newTreeModel.toString());
				
			}
			else if(model.getModelName().equals("strictClockBranchRates")){
				System.out.println(newBranchRateModel.getVariable(0));
				newBranchRateModel = (StrictClockBranchRates) model;
				System.out.println(newBranchRateModel.getVariable(0));
			}
			else if(model.getModelName().equals("SpectrumModel")){
				newSpectrumModel = SpectrumAlignmentModel.duplicateSpectrumAlignmentModel((SpectrumAlignmentModel) model);
			}
			else if(model.getModelName().equals("frequencyModel")){
				f = (FrequencyModel) model;
				System.out.println(Arrays.toString(f.getFrequencies()));
			}
			
		}
		SpectrumAlignmentUtils.compareTwoSpectrumModel(
				spectrumModel, newSpectrumModel);
		SpectrumTreeLikelihood newTreeLikelihood = new SpectrumTreeLikelihood(
				newSpectrumModel, newTreeModel, newSiteModel, null, 
				false, false, true, false, false);
		System.out.println(newTreeLikelihood.getLogLikelihood());
		///////////////////////////////////////////////////////////////////////

		
//		SpectrumAlignmentModel newSpectrumModel = SpectrumAlignmentModel.duplicateSpectrumAlignmentModel(spectrumModel);
		SpectrumTreeLikelihood newSpectrumTreeLikelihood = new SpectrumTreeLikelihood(newSpectrumModel, treeModel,
				newSiteModel, null, false, false, true, false, false);
		double expected = newSpectrumTreeLikelihood.getLogLikelihood();
		double sslikelihood = likelihood.getLogLikelihood();
		assertEquals(expected, sslikelihood, 0);
		
		List<Likelihood> likelihoods = new ArrayList<Likelihood>();
	
		// Prior
		OneOnXPrior oneOnX = new OneOnXPrior();
		oneOnX.addData(popSize);
	
		DistributionLikelihood logNormalLikelihood = new DistributionLikelihood(
				new LogNormalDistribution(1.0, 1.25), 0); 
		logNormalLikelihood.addData(kappa);
	
		likelihoods.add(oneOnX);
		likelihoods.add(logNormalLikelihood);
		likelihoods.add(coalescent);

		Likelihood newPrior = new CompoundLikelihood(0, likelihoods);
		newPrior.setId(CompoundLikelihoodParser.PRIOR);
	
		// Likelihood
		likelihoods.clear();
		likelihoods.add(newTreeLikelihood);
		Likelihood newLikelihood = new CompoundLikelihood(-1, likelihoods);
		newLikelihood.setId(CompoundLikelihoodParser.LIKELIHOOD);
	
		// ShortReadLikelihood
		likelihoods.clear();
		likelihoods.add(newSrpLikelihood);
		Likelihood newShortReadlikelihood = new CompoundLikelihood(-1, likelihoods);
		newShortReadlikelihood.setId(ShortReadsSpectrumLikelihood.SHORT_READ_LIKELIHOOD);
	
		// Posterior
		likelihoods.clear();
		likelihoods.add(newPrior);
		likelihoods.add(newLikelihood);
		likelihoods.add(newSrpLikelihood);
		Likelihood newPosterior = new CompoundLikelihood(0, likelihoods);
		newPosterior.setId(CompoundLikelihoodParser.POSTERIOR);
	
		assertEquals(prior.getLogLikelihood(), newPrior.getLogLikelihood(), 0);
		
		assertEquals(likelihood.getLogLikelihood(), treeLikelihood.getLogLikelihood(), 0);
		assertEquals(newLikelihood.getLogLikelihood(), newTreeLikelihood.getLogLikelihood(), 0);
		
		assertEquals(likelihood.getLogLikelihood(), newLikelihood.getLogLikelihood(), 0);
		assertEquals(treeLikelihood.getLogLikelihood(), newTreeLikelihood.getLogLikelihood(), 0);
		assertEquals(srpLikelihood.getLogLikelihood(), newShortReadlikelihood.getLogLikelihood(), 0);
		assertEquals(posterior.getLogLikelihood(), newPosterior.getLogLikelihood(), 0);
		System.out.println("==============");
	}





}
