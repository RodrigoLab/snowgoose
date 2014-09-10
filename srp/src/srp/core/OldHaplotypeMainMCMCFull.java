package srp.core;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import srp.dr.ext.TreeLikelihoodExt;
//import srp.evolution.haplotypes.old.OldHaplotypeModelUtils;
import srp.evolution.haplotypes.HaplotypeLoggerWithTrueHaplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import srp.evolution.shortreads.ShortReadMapping;
//import srp.likelihood.haplotypes.ShortReadLikelihood;
import dr.evolution.alignment.Alignment;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.tree.TreeLogger;
import dr.evomodel.tree.TreeModel;
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

public class OldHaplotypeMainMCMCFull {

	public static void main(String[] args) throws Exception {

		String dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/testData/";
		int runIndex = 51;
		int totalSamples = 1000;
		int logInterval = 1000;
		int noOfTrueHaplotype = 7;
		int noOfRecoveredHaplotype=7;
		dataDir += "H7_"+runIndex;
		
//		String dataDir = args[0];
//		int runIndex = Integer.parseInt(args[1]);
//		int totalSamples = Integer.parseInt(args[2]);
//		int logInterval = Integer.parseInt(args[3]);
//		int noOfTrueHaplotype = Integer.parseInt(args[4]);
//		int noOfRecoveredHaplotype= Integer.parseInt(args[5]);
		
		String hapRunIndex = "H"+noOfTrueHaplotype+"_"+runIndex;
		String shortReadFile = hapRunIndex +"_Srp.fasta";
		String trueHaplotypeFile = hapRunIndex +"_Srp_fullHaplotype.fasta";
		
		String prefix = dataDir+"FullTree_"+hapRunIndex;
		String logTracerName = prefix+".log";
		String logTreeName = prefix+".trees";
		String logHaplotypeName = prefix+".haplatype";
		String operatorAnalysisFile = prefix+"_operatorAnalysisFile.txt";
		
		DataImporter dataImporter = new DataImporter(dataDir);

//		Alignment shortReads = dataImporter.importShortReads(shortReadFile);
		HaplotypeModel haplotypeModel = new HaplotypeModel(noOfRecoveredHaplotype, 1200);

		Alignment trueAlignment = dataImporter.importAlignment(trueHaplotypeFile);
//		haplotypeModel = new HaplotypeModel(alignmentMapping, trueAlignment);
//		ShortReadLikelihood shortReadLikelihood  = new ShortReadLikelihood(haplotypeModel);
		
		// coalescent
		Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000.0, 100, 100000.0);

		// Random treeModel
		ConstantPopulationModel popModel = new ConstantPopulationModel(popSize, Units.Type.YEARS);
		TreeModel treeModel = MCMCSetupHelperHaplotype.setupRandomTreeModel(popModel, haplotypeModel, Units.Type.YEARS);
		
		// Coalescent likelihood
		CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel,null, new ArrayList<TaxonList>(), popModel);
		coalescent.setId("coalescent");

		// Simulate haplotypes, treeLikelihood
		ShortReadMapping srpMap = null;
		HashMap<String, Object> parameterList = MCMCSetupHelperHaplotype.setupTreeLikelihoodHaplotypeModel(treeModel, haplotypeModel, srpMap);
		Parameter kappa = (Parameter) parameterList.get("kappa");
		Parameter freqs = (Parameter) parameterList.get("freqs");
		StrictClockBranchRates branchRateModel = (StrictClockBranchRates) parameterList.get("branchRateModel");
		TreeLikelihoodExt treeLikelihood = (TreeLikelihoodExt) parameterList.get("treeLikelihood");
		
		// ShortReadLikelihood
//		ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(haplotypeModel);

		// CompoundLikelihood
		HashMap<String, Likelihood> compoundlikelihoods = MCMCSetupHelperHaplotype.setupCompoundLikelihood(
				popSize, kappa, coalescent, treeLikelihood, null);
		Likelihood prior = compoundlikelihoods.get(CompoundLikelihoodParser.PRIOR);
		Likelihood likelihood = compoundlikelihoods.get(CompoundLikelihoodParser.LIKELIHOOD);
//		Likelihood shortReadLikelihood = compoundlikelihoods.get(ShortReadLikelihood.SHORT_READ_LIKELIHOOD);
		Likelihood posterior = compoundlikelihoods.get(CompoundLikelihoodParser.POSTERIOR);
		
		// Operators
		OperatorSchedule schedule = new SimpleOperatorSchedule();
//		ArrayList<MCMCOperator> defalutOperatorsList = 
//		schedule.addOperators(MCMCSetupHelperHaplotype.defalutOperators(haplotypeModel, freqs, popSize, kappa));
//		schedule.addOperators(MCMCSetupHelperHaplotype.defalutTreeOperators(treeModel));
//		MCMCSetupHelperHaplotype.defalutOperatorsOldHaplotype(schedule, haplotypeModel, freqs, popSize, kappa);
		MCMCSetupHelperHaplotype.defalutTreeOperators(schedule, treeModel);
				
		
//		MCMCOperator operator;
//		operator = new RJTreeOperator(haplotypeModel, treeModel);
//		operator.setWeight(100);
//		schedule.addOperator(operator);
		
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
				posterior, prior, likelihood,// shortReadLikelihood,
				rootHeight ,
				//rateParameter,
				popSize, kappa, coalescent
//				freqs
				);
		// System.out
		loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), logInterval, true, logInterval*2);
		MCMCSetupHelperHaplotype.addToLogger(loggers[1],
//				freqs
				posterior, prior, likelihood, //shortReadLikelihood,
				popSize, kappa, coalescent, rootHeight
				);
		// log Tree
		TabDelimitedFormatter treeFormatter = new TabDelimitedFormatter(
				new PrintWriter(new FileOutputStream(new File(logTreeName))));

		loggers[2] = new TreeLogger(treeModel, branchRateModel, null, null,
				treeFormatter, logInterval, true, true, true, null, null);

		// log Haplotype
//		Alignment trueAlignment = dataImporter.importAlignment(trueHaplotypeFile);
//		AlignmentMapping alignmentMapping = new AlignmentMapping(shortReads);
//		ShortReadLikelihood trueSrp = new ShortReadLikelihood(OldHaplotypeModelUtils.factory(shortReads, trueAlignment));
//		System.err.println("\'trueShortReadLikelihood\': "+trueSrp.getLogLikelihood());
		loggers[3] = new HaplotypeLoggerWithTrueHaplotype(haplotypeModel, trueAlignment, logHaplotypeName, logInterval*10);
		
		// MCMC
		MCMCOptions options = MCMCSetupHelperHaplotype.setMCMCOptions(logInterval, totalSamples);
		
		MCMC mcmc = new MCMC("mcmc1");
		mcmc.setShowOperatorAnalysis(true);
		mcmc.setOperatorAnalysisFile(new File(operatorAnalysisFile));
		
		mcmc.init(options, posterior, schedule, loggers);
		mcmc.run();

		
		System.out.println(mcmc.getTimer().toString());
		
	}


}
