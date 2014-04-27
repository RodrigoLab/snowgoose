package test.dr.ext;


import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.logging.LogManager;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.dr.ext.SitePatternsExt;
import srp.dr.ext.TreeLikelihoodExt;
import srp.evolution.haplotypes.HaplotypeModel;
import srp.evolution.shortreads.ShortReadMapping;
import srp.likelihood.haplotypes.ShortReadsHaplotypeLikelihood;
import srp.operator.haplotypes.BaseSingleOperator;
import srp.operator.haplotypes.BasesMultiOperator;
import srp.operator.haplotypes.ColumnOperator;
import srp.operator.haplotypes.HaplotypeRecombinationOperator;
import test.TestUtils;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SitePatterns;
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
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.math.MathUtils;

public class TreeLikelihoodExtTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}


	private GammaSiteModel siteModel;
	private BranchRateModel branchRateModel;
	private Alignment trueAlignment;
	private Tree truePhylogeny;
	
	private ShortReadMapping srpMap;

	@Before
	public void setUp() throws Exception {
//		java.util.logging.ConsoleHandler = NONE;
//		logger.setUseParentHandlers(false)
		
		LogManager.getLogManager().reset();
		
    	// clock model
    	Parameter rateParameter =  new Parameter.Default(StrictClockBranchRates.RATE, 1E-5, 0, 1);
    	branchRateModel = new StrictClockBranchRates(rateParameter);

    	// Sub model
    	Parameter freqs = new Parameter.Default(new double[]{0.25,0.25,0.25,0.25});
    	Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100.0);

    	FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
    	HKY hky = new HKY(kappa, f);

    	//siteModel
    	Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
    	siteModel = new GammaSiteModel(hky);
    	siteModel.setMutationRateParameter(mu);


		String dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/";

		String trueAlignmentFile = "H6_haplotypes.phyml";
		String phylogenyFile = "H6_haplotypes.tree";
		String shortReadFile = "H6_srp.fasta";//"1110_10_align_2.fasta";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		
		trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		truePhylogeny = dataImporter.importTree(phylogenyFile);
		Alignment shortReads = dataImporter.importShortReads(shortReadFile);
		srpMap = new ShortReadMapping(shortReads);
		

	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testUpdatePatternList() throws Exception {
		
		HaplotypeModel haplotypeModel = new HaplotypeModel(trueAlignment);
		
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false);

    	//treeLikelihood
    	SitePatterns patterns = new SitePatterns(trueAlignment, null, 0, -1, 1, true);
    	TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
    			false, false, true, false, false);

    	//treeLikelihoodExt
//    	SitePatternsExt patternsExt = new SitePatternsExt(haplotypeModel, null, 0, -1, 1, true);
        TreeLikelihoodExt treeLikelihoodExt = new TreeLikelihoodExt(haplotypeModel, treeModel, siteModel, branchRateModel, null,
    			false, false, false, false, false);
		// end

        assertEquals(treeLikelihood.getLogLikelihood(), treeLikelihoodExt.getLogLikelihood(), TestUtils.THRESHOLD);
		
        ShortReadsHaplotypeLikelihood srpLikelihood = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);

		SimpleMCMCOperator op = new BasesMultiOperator(haplotypeModel, 10, CoercionMode.COERCION_OFF);

        for (int i = 0; i < 100; i++) {
        	srpLikelihood.storeModelState();
        	op.doOperation();
					
			double likelihood = srpLikelihood.getLogLikelihood();

			
			HaplotypeModel haplotypeModelFull = HaplotypeModel.duplicateHaplotypeModel(haplotypeModel);
			ShortReadsHaplotypeLikelihood likelihoodFull = new ShortReadsHaplotypeLikelihood(haplotypeModelFull, srpMap);
			double fullEvalution =  likelihoodFull.getLogLikelihood();
			
			assertEquals(likelihood, fullEvalution, TestUtils.THRESHOLD);
			
			//treeLikelihood
			treeLikelihoodExt.updatePatternListExt(haplotypeModel);
	    	patterns = new SitePatterns(haplotypeModel, null, 0, -1, 1, true);
	    	treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
	    			false, false, true, false, false);
	    	
    		assertEquals(treeLikelihood.getLogLikelihood(), treeLikelihoodExt.getLogLikelihood(), 1e-8);

    		srpLikelihood.acceptModelState();
        }


	}
	
	@Test
	public void testUpdatePatternList2() throws Exception {
		
		
		String dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/";

		String trueAlignmentFile = "H6_haplotypes.phyml";
		String phylogenyFile = "H6_haplotypes.tree";
		String shortReadFile = "H6_srp.fasta";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		Alignment alignment = dataImporter.importAlignment(trueAlignmentFile);
		
		Tree truePhylogeny = dataImporter.importTree(phylogenyFile);
		
		Alignment shortReads = dataImporter.importShortReads(shortReadFile);
		ShortReadMapping srpMap = new ShortReadMapping(shortReads);
		
		HaplotypeModel trueHaplotypes = new HaplotypeModel(trueAlignment);
		ShortReadsHaplotypeLikelihood srpLikelihood = new ShortReadsHaplotypeLikelihood(trueHaplotypes, srpMap);
	
		
		
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false);
	//		for (int i = 0; i < trueHaplotypes.getHaplotypeCount(); i++) {
	//			treeModel.setTaxonId(i, trueHaplotypes.getTaxonId(i));
	//		}
	//		
	
	//			LikelihoodCalculation likelihoodModel = new LikelihoodCalculation(treeModel, aMap, trueHaplotypes);
	//			likelihoodModel.setPopSize(3000,0,30000);
	////			li.setTreeAndAlignment(treeModel, trueAlignment);
	//			System.out.println(likelihoodModel.getTreeLikelihood());
	//			System.out.println(likelihoodModel.getCoalescentLikelhood());
	//			System.out.println(likelihoodModel.getShortReadsHaplotypeLikelihood());
	//			System.out.println(likelihoodModel.getLoglikelihood());
	//			System.out.println(srpLikelihood.getLogLikelihood());
	//			System.out.println("=== END True Hap===");
				
			
			
	        //CompoundLikelihood
	//	        List<Likelihood> likelihoodsList = new ArrayList<Likelihood>();        
	//	        likelihoodsList.add(coalescent);
	//	        Likelihood prior = new CompoundLikelihood(0, likelihoodsList);
	//	        prior.setId(CompoundLikelihoodParser.PRIOR);
	//
	//	        likelihoodsList.clear();
	//	        likelihoodsList.add(treeLikelihood);
	//	        Likelihood likelihood = new CompoundLikelihood(-1, likelihoodsList);
	//
	//	        likelihoodsList.clear();
	//	        likelihoodsList.add(prior);
	//	        likelihoodsList.add(likelihood);
	//	        Likelihood posterior = new CompoundLikelihood(0, likelihoodsList);
	//	        posterior.setId(CompoundLikelihoodParser.POSTERIOR);
	
				
			// start
		
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

//	    	TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
//	    			false, false, true, false, false);
//	    	treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);

		// end
		
		
		
	//
	//			Haplotypes haplotypes = new Haplotypes(aMap, trueAlignment.getSequenceCount());
    	HaplotypeModel haplotypeModel = trueHaplotypes;
	////			for (int i = 0; i < haplotypes.getHaplotypesCount(); i++) {
	////					haplotypes.randomSeq(i);
	////				System.out.println(haplotypes.getHaplotype(i) );
	////			}
	//			likelihoodModel.updateHaplotypes(haplotypes);
	//			likelihoodModel.calculateShortReadsHaplotypeLikelihoodFull();
	//			System.out.println(likelihoodModel.getTreeLikelihood());
	//			System.out.println(likelihoodModel.getCoalescentLikelhood());
	//			System.out.println(likelihoodModel.getShortReadsHaplotypeLikelihood());
	//			System.out.println(likelihoodModel.getLoglikelihood());
	//			
	//System.out.println("====full====");
	//			 likelihoodModel = new LikelihoodCalculation(treeModel, aMap, haplotypes);
	//				likelihoodModel.setPopSize(3000,0,30000);
	//			likelihoodModel.setPopSize(3000,0,30000);
	////			li.setTreeAndAlignment(treeModel, trueAlignment);
	//			System.out.println(likelihoodModel.getTreeLikelihood());
	//			System.out.println(likelihoodModel.getCoalescentLikelhood());
	//			System.out.println(likelihoodModel.getShortReadsHaplotypeLikelihood());
	//			System.out.println(likelihoodModel.getLoglikelihood());
	//			
	//			
	//
	//			System.out.println("===");
	

		
		
		
		
		srpLikelihood = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		
                
        
        
        int thinning = 1000;
        SitePatternsExt patternsExt = new SitePatternsExt(haplotypeModel, null, 0, -1, 1, true);
        TreeLikelihoodExt treeLikelihoodExt = new TreeLikelihoodExt(haplotypeModel, treeModel, siteModel, branchRateModel, null,
    			false, false, false, false, false);
        
        double likelihood = srpLikelihood.getLogLikelihood() + treeLikelihoodExt.getLogLikelihood();
		likelihood = treeLikelihoodExt.getLogLikelihood();
		System.out.println(likelihood);
		
		BaseSingleOperator op = new BaseSingleOperator(haplotypeModel);
        for (int i = 0; i < 100; i++) {
        	op.doOperation();
//			alignment = haplotypeModel.getAlignment();
			patternsExt.updateAlignment(haplotypeModel);
			treeLikelihoodExt.updatePatternListExt();
		}
		
		likelihood = treeLikelihoodExt.getLogLikelihood();
		System.out.println("updated: "+likelihood);
		treeLikelihoodExt = new TreeLikelihoodExt(haplotypeModel, treeModel, siteModel, branchRateModel, null,
    			false, false, false, false, false);
		likelihood = treeLikelihoodExt.getLogLikelihood();
		System.out.println("updated: "+likelihood);
		
		
	
	//		
	//        for (int i = 0; i < 1e5; i++) {
	//
	//			haplotypeModel.swapBase();
	//			srpLikelihood.updateHaplotypes(haplotypeModel);
	//			
	////				likelihoodModel.updateHaplotypes(haplotypes);
	////				int[] x = haplotypes.getSwapInfo();
	////				if(x[2] != x[3]){
	////					System.out.println(Arrays.toString(haplotypes.getSwapInfo()));
	////					System.out
	////					.println(Arrays.toString(patterns2.getStateFrequencies()));
	////			System.out.println(Arrays.toString(patterns2.getPattern(x[0])));
	////				}
	//			
	//			alignment = haplotypeModel.getAlignment();
	//			patternsExt.updateAlignment(alignment);
	//			treeLikelihoodExt.updatePatternList(patternsExt);
	//			
	//			
	////				if(x[2] != x[3]){
	////				System.out
	////						.println(Arrays.toString(patterns2.getStateFrequencies()));
	////				System.out.println(Arrays.toString(patterns2.getPattern(x[0])));
	////		System.out.println("===");
	////				}
	////				double newL = srpLikelihood.getLogLikelihood() + treeLikelihood2.getLogLikelihood();
	//			double newL = treeLikelihoodExt.getLogLikelihood();
	//
	//			
	//			
	//			boolean accept = criterion.accept(likelihood, newL, hastingsRatio, logr);
	////					System.out.println(likelihood +"\t"+newL +"\t"+ logr[0] +"\t"+  accept);
	//			if(accept){
	//				likelihood = newL;
	////					likelihoodModel.acceptState();
	////					srpLikelihood.acceptState();
	//			}
	//			else{
	//				haplotypeModel.reject();
	////					likelihoodModel.restorreState();
	//				srpLikelihood.restoreState();
	//			}
	//			if (i% thinning == 0){
	//				
	//				String temp = i +"\t"+ likelihood +"\t"+ haplotypeModel.calculateSPS() + 
	//						"\t"+ HaplotypeModel.calculeteSPS(trueHaplotypes, haplotypeModel) + "\n";
	//				System.out.print(temp);
	//			}
	//					
	//		}
	//		
	//		System.out.println("===== Done =====");
	//		System.out.println("Likelihood:\t"+srpLikelihood.getLogLikelihood() 
	//				+"\t"+ haplotypeModel.calculateSPS());
	//		sps = HaplotypeModel.calculeteSPSArray(trueHaplotypes, haplotypeModel);
	//		for (int i = 0; i < sps.length; i++) {
	//				System.out.println(Arrays.toString(sps[i]));
	//		}
	//		}catch (Exception e){
	//			e.printStackTrace();
	//		}
	////				Haplotypes.calculeteSPS(trueHaplotypes, haplotypes);
	//			
	//
	//
	//	}
		}

	@Test
	public void testUpdateHaplotypeModel() throws Exception {
		
		HaplotypeModel haplotypeModel = new HaplotypeModel(trueAlignment);
		
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false);

    	//treeLikelihoodExt
        TreeLikelihoodExt treeLikelihoodExt = new TreeLikelihoodExt(haplotypeModel, treeModel, siteModel, branchRateModel, null,
    			false, false, false, false, false);

        SimpleMCMCOperator op = new BasesMultiOperator(haplotypeModel, 10, CoercionMode.COERCION_OFF);

		System.out.println(treeLikelihoodExt.getLogLikelihood());
		System.out.println(treeLikelihoodExt.time1);
        System.out.println(treeLikelihoodExt.time2);
        
//      
        long timeo = 0;
        long timel = 0;
        double ite = 1e4;
        for (int i = 0; i < ite; i++) {
        	treeLikelihoodExt.storeModelState();
        	long t1 = System.nanoTime();
        	op.doOperation();
        	long t2 = System.nanoTime();
        	treeLikelihoodExt.getLogLikelihood();
        	long t3 = System.nanoTime();
        	treeLikelihoodExt.acceptModelState();
        	timeo+= (t2-t1);
        	timel += (t3-t2);
//        	SitePatterns patterns = new SitePatterns(haplotypeModel, null, 0, -1, 1, true);
//        	TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
//        			false, false, true, false, false);
//    		assertEquals(treeLikelihood.getLogLikelihood(), treeLikelihoodExt.getLogLikelihood(), 1e-8);

        }
        System.out.println(timeo/ite/1e3);
        System.out.println(timel/ite/1e3);

        HaplotypeModel dupHaplotypeModel = HaplotypeModel.duplicateHaplotypeModel(haplotypeModel);
    	SitePatterns patterns = new SitePatterns(dupHaplotypeModel, null, 0, -1, 1, true);
    	TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
    			false, false, true, false, false);
		assertEquals(treeLikelihood.getLogLikelihood(), treeLikelihoodExt.getLogLikelihood(), 1e-8);
	    
	
	}
	

	@Test
	public void testUpdateHaplotypeModelMulitOperators() throws Exception {
		
		HaplotypeModel haplotypeModel = new HaplotypeModel(trueAlignment);
		
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false);

    	//treeLikelihoodExt
        TreeLikelihoodExt treeLikelihoodExt = new TreeLikelihoodExt(haplotypeModel, treeModel, siteModel, branchRateModel, null,
    			false, false, false, false, false);

        OperatorSchedule schedule = new SimpleOperatorSchedule();
		MCMCOperator op; 

		op = new BaseSingleOperator(haplotypeModel);
		schedule.addOperator(op);

		op = new BasesMultiOperator(haplotypeModel, 10, CoercionMode.COERCION_OFF);
		schedule.addOperator(op);

		op = new ColumnOperator(haplotypeModel, haplotypeModel.getHaplotypeCount(), null);
		schedule.addOperator(op);

		op = new HaplotypeRecombinationOperator(haplotypeModel, 10, CoercionMode.COERCION_OFF);
		schedule.addOperator(op);

		double ite = 1e3;
        for (int i = 0; i < ite; i++) {
			
			boolean operatorSucceeded = true;
            boolean accept = false;
            
			final int opIndex = schedule.getNextOperatorIndex();
            final MCMCOperator mcmcOperator = schedule.getOperator(opIndex);
            
            try{
            	mcmcOperator.operate();
            }
			catch (OperatorFailedException e) {
                operatorSucceeded = false;
			}
			
			if(operatorSucceeded){

				SitePatterns patterns = new SitePatterns(haplotypeModel, null, 0, -1, 1, true);
	        	TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
	        			false, false, true, false, false);
	    		assertEquals(treeLikelihood.getLogLikelihood(), treeLikelihoodExt.getLogLikelihood(), 1e-8);
				
				double rand = MathUtils.nextDouble();
				accept = rand>0.5;
			}
			if(accept){
				mcmcOperator.accept(0);
				treeLikelihoodExt.acceptModelState();
			}
			else{
				mcmcOperator.reject();
				treeLikelihoodExt.restoreModelState();
			}

        }

        
        
        ite = 1e5;
        for (int i = 0; i < ite; i++) {
			
			boolean operatorSucceeded = true;
            boolean accept = false;
            
			final int opIndex = schedule.getNextOperatorIndex();
            final MCMCOperator mcmcOperator = schedule.getOperator(opIndex);
            
            try{
            	mcmcOperator.operate();
            }
			catch (OperatorFailedException e) {
                operatorSucceeded = false;
			}
			
			if(operatorSucceeded){
//				double logLikelihoodOperator = treeLikelihoodExt.getLogLikelihood();
				
				double rand = MathUtils.nextDouble();
				accept = rand>0.5;
			}
			if(accept){
				mcmcOperator.accept(0);
				treeLikelihoodExt.acceptModelState();
			}
			else{
				mcmcOperator.reject();
				treeLikelihoodExt.restoreModelState();
			}
        
        }
        HaplotypeModel dupHaplotypeModel = HaplotypeModel.duplicateHaplotypeModel(haplotypeModel);
    	SitePatterns patterns = new SitePatterns(dupHaplotypeModel, null, 0, -1, 1, true);
    	TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
    			false, false, true, false, false);
		assertEquals(treeLikelihood.getLogLikelihood(), treeLikelihoodExt.getLogLikelihood(), 1e-8);
	      
	
	}
}
