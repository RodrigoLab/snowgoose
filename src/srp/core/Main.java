package srp.core;
/*
 -XX:CompileThreshold=50000 -XX:+CITime
-Djava.library.path=/home/sw167/PostdocLarge/Software/BEAST/BEASTv1.7.1/lib -Xms128m -Xmx256m
*/

import java.util.Arrays;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.HaplotypeModelUtils;
import dr.evolution.alignment.Alignment;


public class Main {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {

		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/Stage0/";
//		testReadFiles(dataDir);
//		System.setErr(new PrintStream(new OutputStream() {
//		    @Override
//			public void write(int b) {
//		    }
//		}));		
		dataDir = "/home/sw167/workspace/ABI/";
		dataDir = "/home/sw167/workspace/ABI/data/";
//		String truePhylogenyFile = "H6_005_true_tree.trees";
		String shortReadFile = "H6_srp.fasta";
				
				

		String trueAlignmentFile = "H6_005.fasta";
//		String truePhylogenyFile = "H4_haplotypes.tree";
		String finalAlignmentFile = "finalAlignment";
//		String shortReadFile = "unittest/H6_srp_300.fasta";
		
		DataImporter dataImporter = new DataImporter(dataDir);;
		Alignment finalAlignment = dataImporter.importAlignment(finalAlignmentFile);
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
		
		
		HaplotypeModel tHaplotypeModel = new HaplotypeModel(aMap, trueAlignment);
		HaplotypeModel fHaplotypeModel = new HaplotypeModel(aMap, finalAlignment);
		

		int[][] sps = HaplotypeModelUtils.calculeteSPSArray(tHaplotypeModel, fHaplotypeModel);
		for (int i = 0; i < sps.length; i++) {
				System.out.println(Arrays.toString(sps[i]));
		}
		
		

		
		long time1, time2;
		char[] c = new char[1000000];
		time1 = System.currentTimeMillis();
		for (int i = 0; i < 1e6; i++) {
			c[i] = (char) i;}
		for (int q = 0; q < 100; q++) {
		for (int i = 0; i < 1e6; i++) {
			c[i] = (char) i;
//			char t = c[i];
		}}
		time2 = System.currentTimeMillis();
		System.out.println("\t" + (time2 - time1) + "\t");

		StringBuilder sb = new StringBuilder(1000000);
		for (int i = 0; i < 1000000; i++) {
			sb.append(1);
		}
		time1 = System.currentTimeMillis();
		for (int q = 0; q < 100; q++) {
		for (int i = 0; i < 1e6; i++) {
			sb.setCharAt(i, (char) (i));
//			char t = sb.charAt(i);
		}}
		time2 = System.currentTimeMillis();
		System.out.println("\t" + (time2 - time1) + "\t");


		
		
//		testSim1(args);
		
		
		
		//		testAlignmentMapping();
//		testCompoundLikelihood();

//		testAlignmentMappingLikelihood();
		
	}
//
//	private static void testSim1(String[] args){
////		. 1110_10_org_6.phyml 1110_10_org_6.phyml_phyml_tree.txt 1110_10_align_100.fasta 1		
//		String dataDir = args[0];//"/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/";
//		dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/";
//		
//		String trueAlignmentFile = args[1];//"1110_10_org_6.phyml";
//		String truephylogenyFile = args[2];//"1110_10_org_6.phyml_phyml_tree.txt";
//		String shortReadFile = args[3];//"1110_10_align_100.fasta";//"1110_10_align_2.fasta";
//		
//		truephylogenyFile = "130220_H4_haplotypes.tree";
//		shortReadFile = "130220_H4_srp.fasta";
//		
//		DataImporter dataImporter = new DataImporter(dataDir);
//		Alignment alignment = dataImporter.importAlignment(trueAlignmentFile);
//
//		Tree truePhylogeny = dataImporter.importTree(truephylogenyFile);
//		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);
//		
//		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
//		AlignmentMapping aMap = new AlignmentMapping(shortReads);
//		
//		
//		int[][] sps;
//		final HaplotypeModel trueHaplotypes = new HaplotypeModel(aMap, alignment);
//		HaplotypeModel haplotypeModel;// = new Haplotypes(aMap, alignment);//.getSequenceCount());
////		if(args[4].equals("0")){
////			haplotypes = new Haplotypes(aMap, alignment);//.getSequenceCount());
////		}
////		else{
////			haplotypes = new Haplotypes(aMap, alignment.getSequenceCount());
////		}
//		
////		{Alignment startingAlignment = trueHaplotypes.swapAlignment();
////		haplotypes = new Haplotypes(aMap, startingAlignment);
////		System.out.println(haplotypes.toString());}
//		
//		haplotypeModel = new HaplotypeModel(aMap, alignment);
//		{int max = haplotypeModel.getHaplotypeCount() * haplotypeModel.getHaplotypeLength() *1;
//		System.out.println("Max swap: "+max);
//		for (int i = 0; i < max; i++) {
//			haplotypeModel.swapBase();
//		}}
////		
////		haplotypes = new Haplotypes(aMap, alignment.getSequenceCount());
//		
//	
//		sps = HaplotypeModel.calculeteSPSArray(trueHaplotypes, haplotypeModel);
//		for (int i = 0; i < sps.length; i++) {
//				System.out.println(Arrays.toString(sps[i]));
//		}
//		
//		
//		int ite = Integer.parseInt(args[5]);
//		int thinning = ite/10000;
//		if (thinning < 1000){
//			thinning = 1000;
//		}
//		ite = (int) 1e6;
//		
//
//		
//        //CompoundLikelihood
////        List<Likelihood> likelihoodsList = new ArrayList<Likelihood>();        
////        likelihoodsList.add(coalescent);
////        Likelihood prior = new CompoundLikelihood(0, likelihoodsList);
////        prior.setId(CompoundLikelihoodParser.PRIOR);
////
////        likelihoodsList.clear();
////        likelihoodsList.add(treeLikelihood);
////        Likelihood likelihood = new CompoundLikelihood(-1, likelihoodsList);
////
////        likelihoodsList.clear();
////        likelihoodsList.add(prior);
////        likelihoodsList.add(likelihood);
////        Likelihood posterior = new CompoundLikelihood(0, likelihoodsList);
////        posterior.setId(CompoundLikelihoodParser.POSTERIOR);
//
//			
//		// start
//		try {
//			Path logFileName = Paths.get(dataDir+"likelihood.log");
//			Path resultFileName = Paths.get(dataDir+"result.log");
//			BufferedWriter logFile = Files.newBufferedWriter(logFileName, StandardCharsets.UTF_8);
//			BufferedWriter resultFile = Files.newBufferedWriter(resultFileName, StandardCharsets.UTF_8);
//			logFile.write("ite\tlikelihood\tsps1\tsps2\n");
//			
//			
//			
//	    	Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000,0,10000);
//	    	ConstantPopulationModel constantModel = new ConstantPopulationModel(popSize, Units.Type.DAYS);//createRandomInitialTree(popSize);
//	
//	    	CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel, null, new ArrayList<TaxonList>(), constantModel);
//	    	coalescent.setId("coalescent");
//	
//	    	// clock model
//	    	Parameter rateParameter =  new Parameter.Default(StrictClockBranchRates.RATE, 1, 0, 1);//1e5-5
//	    	StrictClockBranchRates branchRateModel = new StrictClockBranchRates(rateParameter);
//	
//	    	// Sub model
//	    	System.out.println(Arrays.toString(alignment.getStateFrequencies()));
//	    	Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
//	    	Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 2.0, 0, 100.0);
//	
//	    	FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
//	    	HKY hky = new HKY(kappa, f);
//	
//	    	//siteModel
//	    	GammaSiteModel siteModel = new GammaSiteModel(hky);
//	    	Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
//	    	siteModel.setMutationRateParameter(mu);
//	
//	    	SitePatternsExt patterns = new SitePatternsExt(haplotypeModel, null, 0, -1, 1, true);
//
//	    	TreeLikelihoodExt treeLikelihood = new TreeLikelihoodExt(patterns, treeModel, siteModel, branchRateModel, null,
//	  			false, false, true, false, false);
//	    	treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);
//	
//			// end
//
////			SitePatterns patterns2 = new SitePatterns(alignment, null, 0, -1, 1, true);
////	    	TreeLikelihood treeLikelihood2 = new TreeLikelihood(patterns2, treeModel, siteModel, branchRateModel, null,
////	  			false, false, true, false, false);
//			
//			sps = HaplotypeModel.calculeteSPSArray(trueHaplotypes, trueHaplotypes);
//			String tempTrue = "";
//			{ShortReadLikelihood trueSrpLikelihood = new ShortReadLikelihood(aMap, trueHaplotypes);
//			tempTrue += "TrueSrpLikelihood\t" + trueSrpLikelihood.getLogLikelihood() + "\t"
//					+ treeLikelihood.getLogLikelihood() +"\t"+ trueHaplotypes.calculateSPS() +"\n";}
//			for (int i = 0; i < sps.length; i++) {
//				tempTrue+= Arrays.toString(sps[i])+"\n";
//			}
//			tempTrue +="\n"; 
//			resultFile.write(tempTrue);
//			System.out.println(tempTrue);
//			
//	    	
//			
//			MCMCCriterion criterion = new MCMCCriterion();
//	        double hastingsRatio = 0.0;
//	        double[] logr = {-Double.MAX_VALUE};
//
//	        
//	        double likelihood = 0;
//
//	        
//	        
//	        sps = HaplotypeModel.calculeteSPSArray(trueHaplotypes, haplotypeModel);
//			for (int i = 0; i < sps.length; i++) {
//					resultFile.write(Arrays.toString(sps[i])+"\n");
//			}
//
//			
////			alignment = haplotypes.getAlignment();
////			patterns.updateAlignment(alignment);
////			treeLikelihood.updatePatternList(patterns);
////			likelihood =  srpLikelihood.getLogLikelihood();
////			likelihood = treeLikelihood.getLogLikelihood();
//			ShortReadLikelihood srpLikelihood = new ShortReadLikelihood(aMap, haplotypeModel);
//			likelihood =  srpLikelihood.getLogLikelihood() + treeLikelihood.getLogLikelihood();
//			
//			resultFile.write("Likelihood:\t"+likelihood
//					+"\t"+ haplotypeModel.calculateSPS()+"\n");
////			
//	        for (int i = 0; i < ite; i++) {
//        	
//				haplotypeModel.swapBase();
//				srpLikelihood.updateHaplotypes(haplotypeModel);
//
////				alignment = .getAlignment();
//				patterns.updateAlignment(haplotypeModel);
//				treeLikelihood.updatePatternList(patterns);
//				//srpLikelihood.getLogLikelihood();// +
//				
//				double srpL =  srpLikelihood.getLogLikelihood();
//				double treeL = treeLikelihood.getLogLikelihood();
//
//				double newL = srpL + treeL; 
////				double newL = treeL;
//
//				
////				SitePatterns patterns2 = new SitePatterns(alignment, null, 0, -1, 1, true);
////		    	TreeLikelihood treeLikelihood2 = new TreeLikelihood(patterns2, treeModel, siteModel, branchRateModel, null,
////		  			false, false, true, false, false);
////				double checkL = treeLikelihood2.getLogLikelihood();
////				if(newL!=checkL){
////					System.out.println(newL +"\t"+ checkL);
////					System.exit(-1);
////				}
//				
//				boolean accept = criterion.accept(likelihood, newL, hastingsRatio, logr);
//	//				System.out.println(likelihood +"\t"+newL +"\t"+ logr[0] +"\t"+  accept);
//				if(accept){
//					likelihood = newL;
//	//				likelihoodModel.acceptState();
//	//				srpLikelihood.acceptState();
//				}
//				else{
//					haplotypeModel.reject();
//	//				likelihoodModel.restorreState();
//					treeLikelihood.restoreModelState();
//					treeLikelihood.updatePatternList(patterns);
//					
//					srpLikelihood.restoreState();
//				}
//				if (i% thinning == 0){
////					System.out.println(i +"\t"+  likelihood +"\t"+ newL +"\t"+ srpLikelihood.getLogLikelihood());
//					String temp = i +"\t"+ likelihood +"\t"+ haplotypeModel.calculateSPS() + 
//							"\t"+ HaplotypeModel.calculeteSPS(trueHaplotypes, haplotypeModel) + 
//							"\t"+ newL +"\t"+ srpL +"\t"+ treeL +"\n";
//					System.out.print(temp);
//					logFile.write(temp);
//					logFile.flush();
//					
//					resultFile.write("\nIte:\t"+i+"\n");
//					sps = HaplotypeModel.calculeteSPSArray(trueHaplotypes, haplotypeModel);
//					for (int s = 0; s < sps.length; s++) {
//							resultFile.write(Arrays.toString(sps[s])+"\n");
//					}
//					resultFile.flush();
//				}
//						
//			}
//			
//			resultFile.write("\n\n===Done\n\nLikelihood:\t"+likelihood
//					+"\t"+ haplotypeModel.calculateSPS()+"\n");
//			sps = HaplotypeModel.calculeteSPSArray(trueHaplotypes, haplotypeModel);
//			for (int i = 0; i < sps.length; i++) {
//					resultFile.write(Arrays.toString(sps[i])+"\n");
//			}
//			
//			resultFile.write(haplotypeModel.toString()+"\n");
//			
//			logFile.close();
//			resultFile.close();
//		}catch (Exception e){
//			e.printStackTrace();
//		}
////			Haplotypes.calculeteSPS(trueHaplotypes, haplotypes);
//		
////      double likelihood = srpLikelihood.getLogLikelihood() + treeLikelihood.getLogLikelihood();
////		likelihood = treeLikelihood.getLogLikelihood();
////		System.out.println(likelihood);
////		
////		for (int i = 0; i < 100; i++) {
////			haplotypes.swapBase();
////			
////		}
////		alignment = haplotypes.getAlignment();
////		patterns.updateAlignment(alignment);
////		treeLikelihood.updatePatternList(patterns);
//////		
////		likelihood = treeLikelihood.getLogLikelihood();
////		System.out.println("updated: "+likelihood);
////		TreeLikelihood treeLikelihood2 = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
////  			false, false, false, false, false);
////		likelihood = treeLikelihood2.getLogLikelihood();
////		System.out.println("reset  : "+likelihood);
////		System.exit(1);
//	}


	

}
