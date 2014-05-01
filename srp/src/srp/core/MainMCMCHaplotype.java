package srp.core;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import srp.dr.ext.TreeLikelihoodExt;
import srp.evolution.haplotypes.HaplotypeLoggerWithTrueHaplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import srp.evolution.shortreads.ShortReadMapping;
import srp.likelihood.AbstractShortReadsLikelihood;
import srp.likelihood.haplotypes.ShortReadsHaplotypeLikelihood;
import dr.evolution.alignment.Alignment;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.tree.TreeLogger;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inferencexml.model.CompoundLikelihoodParser;

public class MainMCMCHaplotype {

	public static void main(String[] args) throws Exception {

//		String dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/testData/";
//		int runIndex = 51;
//		int totalSamples = 1000;
//		int logInterval = 1000;
//		int noOfTrueHaplotype = 7;
//		int noOfRecoveredHaplotype=7;
//		dataDir += "H7_"+runIndex + File.separator;
//		
//		String dataDir = args[0];
//		int runIndex = Integer.parseInt(args[1]);
//		int totalSamples = Integer.parseInt(args[2]);
//		int logInterval = Integer.parseInt(args[3]);
//		int noOfTrueHaplotype = Integer.parseInt(args[4]);
//		int noOfRecoveredHaplotype= Integer.parseInt(args[5]);
		

		String dataDir;
		int runIndex;
		int totalSamples;
		int logInterval;
		int noOfTrueHaplotype;
		int noOfRecoveredHaplotype;
		boolean randomTree = true;
		boolean randomHaplotype = true;

		
//		boolean commandLine = true;
//		commandLine = false;
		
		if(args.length ==6){
			dataDir = args[0];
			runIndex = Integer.parseInt(args[1]);
			totalSamples = Integer.parseInt(args[2]);
			logInterval = Integer.parseInt(args[3]);
			noOfTrueHaplotype = Integer.parseInt(args[4]);
			noOfRecoveredHaplotype= Integer.parseInt(args[5]);
		}
		
		else{	
			System.out.println("local parameters");
			dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/testData/";
			runIndex = 154;
			dataDir += "H7_"+runIndex+"/";
			
			totalSamples = 200	;
			logInterval = 10000 ;
			
			randomTree = true;
			randomHaplotype = true;

			noOfTrueHaplotype = 7;
			noOfRecoveredHaplotype=7;
		}
		
		String hapRunIndex = "H"+noOfTrueHaplotype+"_"+runIndex;
		String shortReadFile = hapRunIndex +"_Srp.fasta";
		String trueHaplotypeFile = hapRunIndex +"_Srp_fullHaplotype.fasta";
		
		String prefix = dataDir+"FullTree_"+hapRunIndex;
		String logTracerName = prefix+".log";
		String logTreeName = prefix+".trees";
		String logHaplotypeName = prefix+".haplotype";
		String operatorAnalysisFile = prefix+"_operatorAnalysisFile.txt";
		
		
//		
		DataImporter dataImporter = new DataImporter(dataDir);

		Alignment shortReads = dataImporter.importShortReads(shortReadFile);
		ShortReadMapping srpMap = new ShortReadMapping(shortReads);


		
		HaplotypeModel haplotypeModel = null;
		
		if(randomHaplotype){
			haplotypeModel = new HaplotypeModel(noOfRecoveredHaplotype, shortReads.getSiteCount());
		}
		else{
			String partialHaplotypeName = prefix+".haplotypepartial";
			
//			haplotypeModel = dataImporter.importPartialSpectrumFile(partialHaplotypeName);
		}

		// ShortReadLikelihood
		ShortReadsHaplotypeLikelihood srpLikelihood = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);

		
		// coalescent
		Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000.0, 100, 100000.0);

		// Random treeModel
		ConstantPopulationModel popModel = new ConstantPopulationModel(popSize, Units.Type.YEARS);
//		TreeModel treeModel = MCMCSetupHelperHaplotype.setupRandomTreeModel(popModel, haplotypeModel, Units.Type.YEARS);
		TreeModel treeModel;
		if(randomTree){
			treeModel = MCMCSetupHelperHaplotype.setupRandomTreeModel(popModel, haplotypeModel, Units.Type.YEARS);
		}
		else{
			String partialTreeName = prefix+".treespartial";
			Tree partialPhylogeny = dataImporter.importTree(partialTreeName);
			treeModel = new TreeModel(TreeModel.TREE_MODEL, partialPhylogeny, false);
		}
		
		// Coalescent likelihood
		CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel,null, new ArrayList<TaxonList>(), popModel);
		coalescent.setId("coalescent");

		// Simulate haplotypes, treeLikelihood
		HashMap<String, Object> parameterList = MCMCSetupHelperHaplotype.setupTreeLikelihoodHaplotypeModel(treeModel, haplotypeModel);

		Parameter kappa = (Parameter) parameterList.get("kappa");
		Parameter freqs = (Parameter) parameterList.get("freqs");
		StrictClockBranchRates branchRateModel = (StrictClockBranchRates) parameterList.get("branchRateModel");
		TreeLikelihoodExt treeLikelihood = (TreeLikelihoodExt) parameterList.get("treeLikelihood");
				
		// CompoundLikelihood
		HashMap<String, Likelihood> compoundlikelihoods = MCMCSetupHelperHaplotype
				.setupCompoundLikelihood(popSize, kappa, coalescent,
						treeLikelihood, srpLikelihood);
		Likelihood prior = compoundlikelihoods.get(CompoundLikelihoodParser.PRIOR);
		Likelihood likelihood = compoundlikelihoods.get(CompoundLikelihoodParser.LIKELIHOOD);
		Likelihood shortReadLikelihood = compoundlikelihoods.get(AbstractShortReadsLikelihood.SHORT_READ_LIKELIHOOD);
		Likelihood posterior = compoundlikelihoods.get(CompoundLikelihoodParser.POSTERIOR);
		
		// Operators
		OperatorSchedule schedule = new SimpleOperatorSchedule();
		MCMCSetupHelperHaplotype.defalutOperators(schedule, haplotypeModel, freqs, popSize, kappa);
//		MCMCSetupHelperHaplotype.defalutTreeOperators(schedule, treeModel);
				

		
		Parameter rootHeight = treeModel.getRootHeightParameter();
		rootHeight.setId("rootHeight");
		
		double total = 0;
		for (int i = 0; i < schedule.getOperatorCount(); i++) {
			MCMCOperator op= schedule.getOperator(i);
			total += op.getWeight() ;
		}
		System.out.println("totalWeight: "+total);
		

		// MCLogger
		MCLogger[] loggers = new MCLogger[4];
		// log tracer
		loggers[0] = new MCLogger(logTracerName, logInterval, false, 0);
		MCMCSetupHelperHaplotype.addToLogger(loggers[0], 
				posterior, prior, likelihood, 
				shortReadLikelihood,
				rootHeight, 
				//rateParameter,
				popSize, kappa, coalescent,
				freqs
				);
		// System.out
		loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), logInterval, true, logInterval*2);
		MCMCSetupHelperHaplotype.addToLogger(loggers[1],
//				freqs
				posterior, prior, likelihood, 
				shortReadLikelihood,
				popSize, kappa, coalescent, rootHeight
				);
		// log Tree
		TabDelimitedFormatter treeFormatter = new TabDelimitedFormatter(
				new PrintWriter(new FileOutputStream(new File(logTreeName))));

		loggers[2] = new TreeLogger(treeModel, branchRateModel, null, null,
				treeFormatter, logInterval, true, true, true, null, null);

		// log Haplotype
		Alignment trueAlignment = dataImporter.importAlignment(trueHaplotypeFile);
		HaplotypeModel trueHaplotypeModel = new HaplotypeModel(trueAlignment);
		ShortReadsHaplotypeLikelihood trueSrp = new ShortReadsHaplotypeLikelihood(trueHaplotypeModel, srpMap);
		System.out.println(trueSrp.getLogLikelihood());
//		AlignmentMapping alignmentMapping = new AlignmentMapping(shortReads);
//		ShortReadsHaplotypeLikelihood trueSrp = new ShortReadsHaplotypeLikelihood(HaplotypeModel.factory(shortReads, trueAlignment), srpMap);
//		System.err.println("\'trueShortReadLikelihood\': "+trueSrp.getLogLikelihood());
		loggers[3] = new HaplotypeLoggerWithTrueHaplotype(haplotypeModel,
				trueAlignment, logHaplotypeName, (int) (logInterval
						* totalSamples * 0.1));
		// loggers[3] = new TreeLogger(treeModel, branchRateModel, null, null,
//				treeFormatter, logInterval, true, true, true, null, null);
		
		// MCMC
		MCMCOptions options = MCMCSetupHelper.setMCMCOptions(logInterval, totalSamples);
		
		MCMC mcmc = new MCMC("mcmc1");
		mcmc.setShowOperatorAnalysis(true);
		mcmc.setOperatorAnalysisFile(new File(operatorAnalysisFile));
		
		mcmc.init(options, posterior, schedule, loggers);
		mcmc.run();

		
		System.out.println(mcmc.getTimer().toString());
		System.out.println("True Likelihood: "+trueSrp.getLogLikelihood());
		srpLikelihood.makeDirty();
		System.out.println(srpLikelihood.getLogLikelihood());
//		treeModel
		TreeLikelihood reCalTreeLikelihood = new TreeLikelihood(
				haplotypeModel, treeModel, treeLikelihood.getSiteModel(), treeLikelihood.getBranchRateModel(), null,
				false, false, true, false, false);
		System.out.println(reCalTreeLikelihood.getLogLikelihood());
	
	}


}
