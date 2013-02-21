package core;
/*
 -XX:CompileThreshold=50000 -XX:+CITime
 
*/
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import likelihood.LikelihoodCalculation;
import likelihood.ShortReadLikelihood;
import alignment.AlignmentMapping;
import alignment.Haplotypes;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.Sequences;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
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
import dr.ext.SitePatternsExt;
import dr.ext.TreeLikelihoodExt;
import dr.inference.mcmc.MCMCCriterion;
import dr.inference.model.Parameter;


public class Main {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/Stage0/";
//		testReadFiles(dataDir);
//		System.setErr(new PrintStream(new OutputStream() {
//		    @Override
//			public void write(int b) {
//		    }
//		}));		
//		testAlignmentMapping();
//		testCompoundLikelihood();
		testSim1(args);
//		testAlignmentMappingLikelihood();
		
	}

	private static void testSim1(String[] args){
//		. 1110_10_org_6.phyml 1110_10_org_6.phyml_phyml_tree.txt 1110_10_align_100.fasta 1		
		String dataDir = args[0];//"/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/";

		String trueAlignmentFile = args[1];//"1110_10_org_6.phyml";
		String truephylogenyFile = args[2];//"1110_10_org_6.phyml_phyml_tree.txt";
		String shortReadFile = args[3];//"1110_10_align_100.fasta";//"1110_10_align_2.fasta";
		
		
		
		DataImporter dataImporter = new DataImporter(dataDir);
		Alignment alignment = dataImporter.importAlignment(trueAlignmentFile);

		Tree truePhylogeny = dataImporter.importTree(truephylogenyFile);
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);
		
		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
		
		
		int[][] sps;
		final Haplotypes trueHaplotypes = new Haplotypes(aMap, alignment);
		Haplotypes haplotypes;// = new Haplotypes(aMap, alignment);//.getSequenceCount());
//		if(args[4].equals("0")){
//			haplotypes = new Haplotypes(aMap, alignment);//.getSequenceCount());
//		}
//		else{
//			haplotypes = new Haplotypes(aMap, alignment.getSequenceCount());
//		}
		
//		{Alignment startingAlignment = trueHaplotypes.swapAlignment();
//		haplotypes = new Haplotypes(aMap, startingAlignment);
//		System.out.println(haplotypes.toString());}
		
		haplotypes = new Haplotypes(aMap, alignment);
		{int max = haplotypes.getHaplotypesCount() * haplotypes.getLength() *1;
		System.out.println("Max swap: "+max);
		for (int i = 0; i < max; i++) {
			haplotypes.swapBase();
		}}
//		
//		haplotypes = new Haplotypes(aMap, alignment.getSequenceCount());
		
	
		sps = Haplotypes.calculeteSPSArray(trueHaplotypes, haplotypes);
		for (int i = 0; i < sps.length; i++) {
				System.out.println(Arrays.toString(sps[i]));
		}
		
		
		int ite = Integer.parseInt(args[5]);
		int thinning = ite/10000;
		if (thinning < 1000){
			thinning = 1000;
		}
		ite = (int) 1e6;
		

		
        //CompoundLikelihood
//        List<Likelihood> likelihoodsList = new ArrayList<Likelihood>();        
//        likelihoodsList.add(coalescent);
//        Likelihood prior = new CompoundLikelihood(0, likelihoodsList);
//        prior.setId(CompoundLikelihoodParser.PRIOR);
//
//        likelihoodsList.clear();
//        likelihoodsList.add(treeLikelihood);
//        Likelihood likelihood = new CompoundLikelihood(-1, likelihoodsList);
//
//        likelihoodsList.clear();
//        likelihoodsList.add(prior);
//        likelihoodsList.add(likelihood);
//        Likelihood posterior = new CompoundLikelihood(0, likelihoodsList);
//        posterior.setId(CompoundLikelihoodParser.POSTERIOR);

			
		// start
		try {
			Path logFileName = Paths.get(dataDir+"likelihood.log");
			Path resultFileName = Paths.get(dataDir+"result.log");
			BufferedWriter logFile = Files.newBufferedWriter(logFileName, StandardCharsets.UTF_8);
			BufferedWriter resultFile = Files.newBufferedWriter(resultFileName, StandardCharsets.UTF_8);
			logFile.write("ite\tlikelihood\tsps1\tsps2\n");
			
			
			
	    	Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000,0,10000);
	    	ConstantPopulationModel constantModel = new ConstantPopulationModel(popSize, Units.Type.DAYS);//createRandomInitialTree(popSize);
	
	    	CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel, null, new ArrayList<TaxonList>(), constantModel);
	    	coalescent.setId("coalescent");
	
	    	// clock model
	    	Parameter rateParameter =  new Parameter.Default(StrictClockBranchRates.RATE, 1, 0, 1);//1e5-5
	    	StrictClockBranchRates branchRateModel = new StrictClockBranchRates(rateParameter);
	
	    	// Sub model
	    	System.out.println(Arrays.toString(alignment.getStateFrequencies()));
	    	Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
	    	Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 2.0, 0, 100.0);
	
	    	FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
	    	HKY hky = new HKY(kappa, f);
	
	    	//siteModel
	    	GammaSiteModel siteModel = new GammaSiteModel(hky);
	    	Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
	    	siteModel.setMutationRateParameter(mu);
	
	    	SitePatternsExt patterns = new SitePatternsExt(haplotypes, null, 0, -1, 1, true);

	    	TreeLikelihoodExt treeLikelihood = new TreeLikelihoodExt(patterns, treeModel, siteModel, branchRateModel, null,
	  			false, false, true, false, false);
	    	treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);
	
			// end

//			SitePatterns patterns2 = new SitePatterns(alignment, null, 0, -1, 1, true);
//	    	TreeLikelihood treeLikelihood2 = new TreeLikelihood(patterns2, treeModel, siteModel, branchRateModel, null,
//	  			false, false, true, false, false);
			
			sps = Haplotypes.calculeteSPSArray(trueHaplotypes, trueHaplotypes);
			String tempTrue = "";
			{ShortReadLikelihood trueSrpLikelihood = new ShortReadLikelihood(aMap, trueHaplotypes);
			tempTrue += "TrueSrpLikelihood\t" + trueSrpLikelihood.getLogLikelihood() + "\t"
					+ treeLikelihood.getLogLikelihood() +"\t"+ trueHaplotypes.calculateSPS() +"\n";}
			for (int i = 0; i < sps.length; i++) {
				tempTrue+= Arrays.toString(sps[i])+"\n";
			}
			tempTrue +="\n"; 
			resultFile.write(tempTrue);
			System.out.println(tempTrue);
			
	    	
			
			MCMCCriterion criterion = new MCMCCriterion();
	        double hastingsRatio = 0.0;
	        double[] logr = {-Double.MAX_VALUE};

	        
	        double likelihood = 0;

	        
	        
	        sps = Haplotypes.calculeteSPSArray(trueHaplotypes, haplotypes);
			for (int i = 0; i < sps.length; i++) {
					resultFile.write(Arrays.toString(sps[i])+"\n");
			}

			
//			alignment = haplotypes.getAlignment();
//			patterns.updateAlignment(alignment);
//			treeLikelihood.updatePatternList(patterns);
//			likelihood =  srpLikelihood.getLogLikelihood();
//			likelihood = treeLikelihood.getLogLikelihood();
			ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(aMap, haplotypes);
			likelihood =  srpLikelihood.getLogLikelihood() + treeLikelihood.getLogLikelihood();
			
			resultFile.write("Likelihood:\t"+likelihood
					+"\t"+ haplotypes.calculateSPS()+"\n");
//			
	        for (int i = 0; i < ite; i++) {
        	
				haplotypes.swapBase();
				srpLikelihood.updateHaplotypes(haplotypes);

				alignment = haplotypes.getAlignment();
				patterns.updateAlignment(alignment);
				treeLikelihood.updatePatternList(patterns);
				//srpLikelihood.getLogLikelihood();// +
				
				double srpL =  srpLikelihood.getLogLikelihood();
				double treeL = treeLikelihood.getLogLikelihood();

				double newL = srpL + treeL; 
//				double newL = treeL;

				
//				SitePatterns patterns2 = new SitePatterns(alignment, null, 0, -1, 1, true);
//		    	TreeLikelihood treeLikelihood2 = new TreeLikelihood(patterns2, treeModel, siteModel, branchRateModel, null,
//		  			false, false, true, false, false);
//				double checkL = treeLikelihood2.getLogLikelihood();
//				if(newL!=checkL){
//					System.out.println(newL +"\t"+ checkL);
//					System.exit(-1);
//				}
				
				boolean accept = criterion.accept(likelihood, newL, hastingsRatio, logr);
	//				System.out.println(likelihood +"\t"+newL +"\t"+ logr[0] +"\t"+  accept);
				if(accept){
					likelihood = newL;
	//				likelihoodModel.acceptState();
	//				srpLikelihood.acceptState();
				}
				else{
					haplotypes.reject();
	//				likelihoodModel.restorreState();
					treeLikelihood.restoreModelState();
					treeLikelihood.updatePatternList(patterns);
					
					srpLikelihood.restoreState();
				}
				if (i% thinning == 0){
//					System.out.println(i +"\t"+  likelihood +"\t"+ newL +"\t"+ srpLikelihood.getLogLikelihood());
					String temp = i +"\t"+ likelihood +"\t"+ haplotypes.calculateSPS() + 
							"\t"+ Haplotypes.calculeteSPS(trueHaplotypes, haplotypes) + 
							"\t"+ newL +"\t"+ srpL +"\t"+ treeL +"\n";
					System.out.print(temp);
					logFile.write(temp);
					logFile.flush();
					
					resultFile.write("\nIte:\t"+i+"\n");
					sps = Haplotypes.calculeteSPSArray(trueHaplotypes, haplotypes);
					for (int s = 0; s < sps.length; s++) {
							resultFile.write(Arrays.toString(sps[s])+"\n");
					}
					resultFile.flush();
				}
						
			}
			
			resultFile.write("\n\n===Done\n\nLikelihood:\t"+likelihood
					+"\t"+ haplotypes.calculateSPS()+"\n");
			sps = Haplotypes.calculeteSPSArray(trueHaplotypes, haplotypes);
			for (int i = 0; i < sps.length; i++) {
					resultFile.write(Arrays.toString(sps[i])+"\n");
			}
			
			resultFile.write(haplotypes.toString()+"\n");
			
			logFile.close();
			resultFile.close();
		}catch (Exception e){
			e.printStackTrace();
		}
//			Haplotypes.calculeteSPS(trueHaplotypes, haplotypes);
		
//      double likelihood = srpLikelihood.getLogLikelihood() + treeLikelihood.getLogLikelihood();
//		likelihood = treeLikelihood.getLogLikelihood();
//		System.out.println(likelihood);
//		
//		for (int i = 0; i < 100; i++) {
//			haplotypes.swapBase();
//			
//		}
//		alignment = haplotypes.getAlignment();
//		patterns.updateAlignment(alignment);
//		treeLikelihood.updatePatternList(patterns);
////		
//		likelihood = treeLikelihood.getLogLikelihood();
//		System.out.println("updated: "+likelihood);
//		TreeLikelihood treeLikelihood2 = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
//  			false, false, false, false, false);
//		likelihood = treeLikelihood2.getLogLikelihood();
//		System.out.println("reset  : "+likelihood);
//		System.exit(1);
	}


	

	private static void testCompoundLikelihood(){
		
		
		
		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/";

		String trueAlignmentFile = "1110_10_org_6.phyml";
		String phylogenyFile = "1110_10_org_6.phyml_phyml_tree.txt";
		String shortReadFile = "1110_10_align_100.fasta";//"1110_10_align_2.fasta";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);

		Tree truePhylogeny = dataImporter.importTree(phylogenyFile);
		
		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
		
		Haplotypes trueHaplotypes = new Haplotypes(aMap, trueAlignment);
		ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(aMap, trueHaplotypes);
	
		Alignment alignment = trueAlignment;


//		
//		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);
//		LikelihoodCalculation li = new LikelihoodCalculation(treeModel, trueAlignment, shortReads);
//		li.setPopSize(3000,0,30000);
//		li.setTreeAndAlignment(treeModel, trueAlignment);

		
//		System.out.println("Likelihood:\t"+srpLikelihood.getLogLikelihood() 
//				+"\t"+ trueHaplotypes.calculateSPS());
		int[][] sps = Haplotypes.calculeteSPSArray(trueHaplotypes, trueHaplotypes);
//		for (int i = 0; i < sps.length; i++) {
//				System.out.println(Arrays.toString(sps[i]));
//		}
		
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);


//		LikelihoodCalculation likelihoodModel = new LikelihoodCalculation(treeModel, aMap, trueHaplotypes);
//		likelihoodModel.setPopSize(3000,0,30000);
////		li.setTreeAndAlignment(treeModel, trueAlignment);
//		System.out.println(likelihoodModel.getTreeLikelihood());
//		System.out.println(likelihoodModel.getCoalescentLikelhood());
//		System.out.println(likelihoodModel.getShortReadLikelihood());
//		System.out.println(likelihoodModel.getLoglikelihood());
//		System.out.println(srpLikelihood.getLogLikelihood());
//		System.out.println("=== END True Hap===");
			
		
		
        //CompoundLikelihood
//        List<Likelihood> likelihoodsList = new ArrayList<Likelihood>();        
//        likelihoodsList.add(coalescent);
//        Likelihood prior = new CompoundLikelihood(0, likelihoodsList);
//        prior.setId(CompoundLikelihoodParser.PRIOR);
//
//        likelihoodsList.clear();
//        likelihoodsList.add(treeLikelihood);
//        Likelihood likelihood = new CompoundLikelihood(-1, likelihoodsList);
//
//        likelihoodsList.clear();
//        likelihoodsList.add(prior);
//        likelihoodsList.add(likelihood);
//        Likelihood posterior = new CompoundLikelihood(0, likelihoodsList);
//        posterior.setId(CompoundLikelihoodParser.POSTERIOR);

			
		// start
		try{
    	Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000,0,10000);
    	ConstantPopulationModel constantModel = new ConstantPopulationModel(popSize, Units.Type.DAYS);//createRandomInitialTree(popSize);

    	CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel, null, new ArrayList<TaxonList>(), constantModel);
    	coalescent.setId("coalescent");

    	// clock model
    	Parameter rateParameter =  new Parameter.Default(StrictClockBranchRates.RATE, 1E-5, 0, 1);
    	StrictClockBranchRates branchRateModel = new StrictClockBranchRates(rateParameter);

    	// Sub model
    	Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
    	Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100.0);

    	FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
    	HKY hky = new HKY(kappa, f);

    	//siteModel
    	GammaSiteModel siteModel = new GammaSiteModel(hky);
    	Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
    	siteModel.setMutationRateParameter(mu);

    	//treeLikelihood
    	SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

//    	TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
//    			false, false, true, false, false);
//    	treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);

		// end
		
		
		
//
//		Haplotypes haplotypes = new Haplotypes(aMap, trueAlignment.getSequenceCount());
    	Haplotypes haplotypes = trueHaplotypes;
////		for (int i = 0; i < haplotypes.getHaplotypesCount(); i++) {
////				haplotypes.randomSeq(i);
////			System.out.println(haplotypes.getHaplotype(i) );
////		}
//		likelihoodModel.updateHaplotypes(haplotypes);
//		likelihoodModel.calculateShortReadLikelihoodFull();
//		System.out.println(likelihoodModel.getTreeLikelihood());
//		System.out.println(likelihoodModel.getCoalescentLikelhood());
//		System.out.println(likelihoodModel.getShortReadLikelihood());
//		System.out.println(likelihoodModel.getLoglikelihood());
//		
//System.out.println("====full====");
//		 likelihoodModel = new LikelihoodCalculation(treeModel, aMap, haplotypes);
//			likelihoodModel.setPopSize(3000,0,30000);
//		likelihoodModel.setPopSize(3000,0,30000);
////		li.setTreeAndAlignment(treeModel, trueAlignment);
//		System.out.println(likelihoodModel.getTreeLikelihood());
//		System.out.println(likelihoodModel.getCoalescentLikelhood());
//		System.out.println(likelihoodModel.getShortReadLikelihood());
//		System.out.println(likelihoodModel.getLoglikelihood());
//		
//		
//
//		System.out.println("===");


		
		System.out.println("Likelihood:\t"+srpLikelihood.getLogLikelihood() 
				+"\t"+ haplotypes.calculateSPS());
		sps = Haplotypes.calculeteSPSArray(trueHaplotypes, haplotypes);
		for (int i = 0; i < sps.length; i++) {
				System.out.println(Arrays.toString(sps[i]));
		
		}
		
		
		srpLikelihood = new ShortReadLikelihood(aMap, haplotypes);
		
        
		MCMCCriterion criterion = new MCMCCriterion();
        double hastingsRatio = 0.0;
        double[] logr = {-Double.MAX_VALUE};
        
        
        
        
        
        
        int thinning = 1000;
        SitePatternsExt patterns2 = new SitePatternsExt(alignment, null, 0, -1, 1, true);
//        final TreeLikelihood treeLikelihood2 = new TreeLikelihood(patterns2, treeModel, siteModel, branchRateModel, null,
//    			false, false, true, false, false);
//	
        TreeLikelihoodExt treeLikelihood2 = new TreeLikelihoodExt(patterns2, treeModel, siteModel, branchRateModel, null,
    			false, false, false, false, false);
        double likelihood = srpLikelihood.getLogLikelihood() + treeLikelihood2.getLogLikelihood();
		likelihood = treeLikelihood2.getLogLikelihood();
		System.out.println(likelihood);
		
		for (int i = 0; i < 100; i++) {
			haplotypes.swapBase();
			alignment = haplotypes.getAlignment();
			patterns2.updateAlignment(alignment);
			treeLikelihood2.updatePatternList(patterns2);
		}
		
		likelihood = treeLikelihood2.getLogLikelihood();
		System.out.println("updated: "+likelihood);
		treeLikelihood2 = new TreeLikelihoodExt(patterns2, treeModel, siteModel, branchRateModel, null,
    			false, false, false, false, false);
		likelihood = treeLikelihood2.getLogLikelihood();
		System.out.println("updated: "+likelihood);
//		System.exit(1);
        for (int i = 0; i < 1e5; i++) {

			haplotypes.swapBase();
			srpLikelihood.updateHaplotypes(haplotypes);
			
//			likelihoodModel.updateHaplotypes(haplotypes);
//			int[] x = haplotypes.getSwapInfo();
//			if(x[2] != x[3]){
//				System.out.println(Arrays.toString(haplotypes.getSwapInfo()));
//				System.out
//				.println(Arrays.toString(patterns2.getStateFrequencies()));
//		System.out.println(Arrays.toString(patterns2.getPattern(x[0])));
//			}
			
			alignment = haplotypes.getAlignment();
			patterns2.updateAlignment(alignment);
			treeLikelihood2.updatePatternList(patterns2);
			
			
//			if(x[2] != x[3]){
//			System.out
//					.println(Arrays.toString(patterns2.getStateFrequencies()));
//			System.out.println(Arrays.toString(patterns2.getPattern(x[0])));
//	System.out.println("===");
//			}
//			double newL = srpLikelihood.getLogLikelihood() + treeLikelihood2.getLogLikelihood();
			double newL = treeLikelihood2.getLogLikelihood();

			
			
			boolean accept = criterion.accept(likelihood, newL, hastingsRatio, logr);
//				System.out.println(likelihood +"\t"+newL +"\t"+ logr[0] +"\t"+  accept);
			if(accept){
				likelihood = newL;
//				likelihoodModel.acceptState();
//				srpLikelihood.acceptState();
			}
			else{
				haplotypes.reject();
//				likelihoodModel.restorreState();
				srpLikelihood.restoreState();
			}
			if (i% thinning == 0){
				
				String temp = i +"\t"+ likelihood +"\t"+ haplotypes.calculateSPS() + 
						"\t"+ Haplotypes.calculeteSPS(trueHaplotypes, haplotypes) + "\n";
				System.out.print(temp);
			}
					
		}
		
		System.out.println("===== Done =====");
		System.out.println("Likelihood:\t"+srpLikelihood.getLogLikelihood() 
				+"\t"+ haplotypes.calculateSPS());
		sps = Haplotypes.calculeteSPSArray(trueHaplotypes, haplotypes);
		for (int i = 0; i < sps.length; i++) {
				System.out.println(Arrays.toString(sps[i]));
		}
		}catch (Exception e){
			e.printStackTrace();
		}
//			Haplotypes.calculeteSPS(trueHaplotypes, haplotypes);
		

	}



	private static void testAlignmentMappingLikelihood(){
		
		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/";

		String trueAlignmentFile = "1110_10_org.phyml";
//		String trueAlignmentFile = "zz.fasta";
		String phylogenyFile = "1110_10_org.phyml_phyml_tree.txt";
		String shortReadFile = "1110_10_align_100.fasta";//"1110_10_align_2.fasta";
		String refSeqFile = "1110_10.ref";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		Tree truePhylogeny = dataImporter.importTree(phylogenyFile);
//		Sequences shortReads = dataImporter.importSequence(shortReadFile);
		Sequence refSeq = dataImporter.importRefSeq(refSeqFile);
		
		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
	
		Path logFileName = Paths.get(dataDir+"likelihood.log");
		try {
			    
			BufferedWriter logFile = Files.newBufferedWriter(logFileName, StandardCharsets.UTF_8);
	//		AlignmentMatrix haplotypes = new AlignmentMatrix(aMap, 20);
			Haplotypes haplotypes = new Haplotypes(aMap, trueAlignment);
			ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(aMap, haplotypes);
			System.out.println("True Alignment Likelihood:\t"+srpLikelihood.getLogLikelihood() +"\t"+ haplotypes.calculateSPS());
	
			haplotypes = new Haplotypes(aMap, trueAlignment.getSequenceCount());
			for (int i = 0; i < haplotypes.getHaplotypesCount(); i++) {
				haplotypes.randomSeq(i);
//				System.out.println(haplotypes.getHaplotype(i) );
			}
	
			srpLikelihood = new ShortReadLikelihood(aMap, haplotypes);
			System.out.println("Likelihood:\t"+srpLikelihood.getLogLikelihood() +"\t"+ haplotypes.calculateSPS());
			
			double likelihood = srpLikelihood.getLogLikelihood(); 
			
			MCMCCriterion criterion = new MCMCCriterion();
	        double hastingsRatio = 0.0;
	        double[] logr = {-Double.MAX_VALUE};
	        
			for (int i = 0; i < 1e5; i++) {
	
				haplotypes.swapBase();
	
				srpLikelihood.updateHaplotypes(haplotypes);
				
	//			ShortReadLikelihood srL_temp = new ShortReadLikelihood(aMap, haplotypes);
				double newL = srpLikelihood.getLogLikelihood();
				
				boolean accept = criterion.accept(likelihood, newL, hastingsRatio, logr);
//				System.out.println(likelihood +"\t"+newL +"\t"+ logr[0] +"\t"+  accept);
				if(accept){
					likelihood = newL;
				}
				else{
					haplotypes.reject();
					
				}
				double sps = haplotypes.calculateSPS();
				logFile.write(i +"\t"+ likelihood +"\t"+ sps +"\n");
					
			}
			List<String> al = new ArrayList<String>();

			al.add(srpLikelihood.getLogLikelihood()+"\n");
			al.add(haplotypes.toString());
			String s = "";
			Files.write(Paths.get(dataDir+"tempAlignment"), al, StandardCharsets.UTF_8);
			
//			System.out.println(haplotypes.toString());
//			srpLikelihood = new ShortReadLikelihood(aMap, haplotypes);
			System.out.println("Likelihood:\t"+srpLikelihood.getLogLikelihood() +"\t"+ haplotypes.calculateSPS());
			

		} catch (IOException x) {
		    System.err.format("IOException: %s%n", x);
		}
	}

	private static void testAlignmentMapping() {
		String dir = "/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/";
		String refFileName = "1110_10.ref";
		String srpFileName = "1110_10_align_test.fasta";
		
		try {
			

			DataImporter dataImporter = new DataImporter(dir);
//			Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
//			Tree truePhylogeny = dataImporter.importTree(phylogenyFile);
			Alignment shortReads = dataImporter.importAlignment(srpFileName);
			Sequence refSeq = dataImporter.importRefSeq(refFileName);
			System.out.println("===");
			System.out.println(shortReads.getPatternCount());
			System.out.println(shortReads.getSiteCount());
			System.out.println(shortReads.getStateCount());
			System.out.println(shortReads.getTaxonCount());
//			System.out.println(shortReads.get
			for (int i = 0; i < shortReads.getSequenceCount(); i++) {
				Sequence s = shortReads.getSequence(i);
				Taxon t = s.getTaxon();
//				System.out.println(s.getId() +"\t"+ t.getId());	
//				System.out.println(s.getChar(0));
//				System.out.println(s.getLength());
				
				
				
//				System.out.println(s.getLength() +"\t"+ s.getSequenceString());
			}
//			BufferedReader in = new BufferedReader(new FileReader(dir+refFileName));
//			in.readLine();
////			String refSeq = in.readLine();
//			
//			in = new BufferedReader(new FileReader(dir+srpFileName));
//			String input = in.readLine();
//			in.readLine();
//			int count = 0;
//			
//			LinkedHashMap<String, String> alignment = new LinkedHashMap();
//			while((input=in.readLine())!= null){
//				if(input.indexOf(">")!= -1){
//					alignment.put(input.substring(1), in.readLine());
//				}
//				count++;
//				if(count==9){
//					break;
//				}
//			}
//			
//			System.out.println(alignment.toString() );
			
			AlignmentMapping amap = new AlignmentMapping(shortReads);
//			AlignmentMapping amap = new AlignmentMapping(1200);
//			for (int i = 0; i < shortReads.getSequenceCount(); i++) {
//				Sequence s = shortReads.getSequence(i);
//				amap.addSequence(s);
//			}
//			shortReads.
//			
//			for (String key: alignment.keySet()){
////				System.out.println(key);
//				
////				System.out.println(index);// +"\t"+ alignment.get(key));
//				amap.addSeq(key, alignment.get(key));
//			}
//			System.out.println(am.toString());
			Haplotypes haplotypes = new Haplotypes(amap, 5);
//			alignmentMatrix.testGetSeq();
//			
//			alignmentMatrix.testMultipleSeq();
			haplotypes.getAlignment();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}

	private static void testReadFiles(String dataDir) {
		
		String trueAlignmentFile = "121101_true_seqs.fasta";
//		String trueAlignmentFile = "zz.fasta";
		String truePhylogenyFile = "121101_true_tree.newick";
		String shortReadFile = "121101_short_reads_10.fasta";
		String refSeqFile = "121101_ref.fasta";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		Tree truePhylogeny = dataImporter.importTree(truePhylogenyFile);
		Sequences shortReads = dataImporter.importSequence(shortReadFile);
		Sequence refSeq = dataImporter.importRefSeq(refSeqFile);
		
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);
		
		
		LikelihoodCalculation li = new LikelihoodCalculation(treeModel, trueAlignment, shortReads);
		li.setPopSize(3000,0,30000);
//		li.setTreeAndAlignment(treeModel, trueAlignment);
		System.out.println(li.getTreeLikelihood());
		System.out.println(li.getCoalescentLikelhood());

		System.out.println(li.getShortReadLikelihood());
		
	}

}
