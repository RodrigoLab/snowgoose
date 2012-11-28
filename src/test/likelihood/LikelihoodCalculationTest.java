package test.likelihood;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import likelihood.LikelihoodCalculation;
import likelihood.ShortReadLikelihood;

import operator.AlignmentOperator;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

//import test.dr.integration.PathSampling;

import core.DataImporter;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.Sequences;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.tree.TreeModelParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.Scalable;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.operators.UniformOperator;
import dr.inference.operators.UpDownOperator;
import dr.inference.trace.ArrayTraceList;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceCorrelation;
import dr.inferencexml.model.CompoundLikelihoodParser;

public class LikelihoodCalculationTest {

	private static final String[][] PRIMATES_TAXON_SEQUENCE = {{"human", "chimp", "bonobo", "gorilla", "orangutan", "siamang"},
	{"AGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGGAGCTTAAACCCCCTTATTTCTACTAGGACTATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAAGAAATTTAGGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTG-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGACCAATGGGACTTAAACCCACAAACACTTAGTTAACAGCTAAGCACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCGGAGCTTGGTAAAAAGAGGCCTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGGCCTCCATGACTTTTTCAAAAGGTATTAGAAAAACCATTTCATAACTTTGTCAAAGTTAAATTATAGGCT-AAATCCTATATATCTTA-CACTGTAAAGCTAACTTAGCATTAACCTTTTAAGTTAAAGATTAAGAGAACCAACACCTCTTTACAGTGA",
	 "AGAAATATGTCTGATAAAAGAATTACTTTGATAGAGTAAATAATAGGAGTTCAAATCCCCTTATTTCTACTAGGACTATAAGAATCGAACTCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTACACCCTTCCCGTACTAAGAAATTTAGGTTAAGCACAGACCAAGAGCCTTCAAAGCCCTCAGCAAGTTA-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATTAATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCAGAGCTTGGTAAAAAGAGGCTTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCTAAAGCTGGTTTCAAGCCAACCCCATGACCTCCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATTACAGGTT-AACCCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGATTAAGAGGACCGACACCTCTTTACAGTGA",
	 "AGAAATATGTCTGATAAAAGAATTACTTTGATAGAGTAAATAATAGGAGTTTAAATCCCCTTATTTCTACTAGGACTATGAGAGTCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAAGAAATTTAGGTTAAACACAGACCAAGAGCCTTCAAAGCTCTCAGTAAGTTA-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATTAATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAATCAGC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTTGAATTTGCAATTCAATATGAAAA-TCACCTCAGAGCTTGGTAAAAAGAGGCTTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCTAAAGCTGGTTTCAAGCCAACCCCATGACCCCCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATTACAGGTT-AAACCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGATTAAGAGGACCAACACCTCTTTACAGTGA",
	 "AGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGAGGTTTAAACCCCCTTATTTCTACTAGGACTATGAGAATTGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTGTCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTCACATCCTTCCCGTACTAAGAAATTTAGGTTAAACATAGACCAAGAGCCTTCAAAGCCCTTAGTAAGTTA-CAACACTTAATTTCTGTAAGGACTGCAAAACCCTACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATCAATGGGACTCAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAGTCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAT-TCACCTCGGAGCTTGGTAAAAAGAGGCCCAGCCTCTGTCTTTAGATTTACAGTCCAATGCCTTA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGACCTTCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAGGTTAAATTACGGGTT-AAACCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCGTTAACCTTTTAAGTTAAAGATTAAGAGTATCGGCACCTCTTTGCAGTGA",
	 "AGAAATATGTCTGACAAAAGAGTTACTTTGATAGAGTAAAAAATAGAGGTCTAAATCCCCTTATTTCTACTAGGACTATGGGAATTGAACCCACCCCTGAGAATCCAAAATTCTCCGTGCCACCCATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTACACCCTTCCCGTACTAAGAAATTTAGGTTA--CACAGACCAAGAGCCTTCAAAGCCCTCAGCAAGTCA-CAGCACTTAATTTCTGTAAGGACTGCAAAACCCCACTTTGCATCAACTGAGCGCAAATCAGCCACTTTAATTAAGCTAAGCCCTCCTAGACCGATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAGTCAAT-TGGCTTCAGTCCAAAGCCCCGGCAGGCCTTAAAGCTGCTCCTTCGAATTTGCAATTCAACATGACAA-TCACCTCAGGGCTTGGTAAAAAGAGGTCTGACCCCTGTTCTTAGATTTACAGCCTAATGCCTTAACTCGGCCATTTTACCGCAAAAAAGGAAGGAATCGAACCTCCTAAAGCTGGTTTCAAGCCAACCCCATAACCCCCATGACTTTTTCAAAAGGTACTAGAAAAACCATTTCGTAACTTTGTCAAAGTTAAATTACAGGTC-AGACCCTGTGTATCTTA-CATTGCAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGACTAAGAGAACCAGCCTCTCTTTGCAATGA",
	 "AGAAATACGTCTGACGAAAGAGTTACTTTGATAGAGTAAATAACAGGGGTTTAAATCCCCTTATTTCTACTAGAACCATAGGAGTCGAACCCATCCTTGAGAATCCAAAACTCTCCGTGCCACCCGTCGCACCCTGTTCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCATACTAAGAAATTTAGGTTAAACACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTAACAAAACTTAATTTCTGCAAGGGCTGCAAAACCCTACTTTGCATCAACCGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATCGATGGGACTTAAACCCATAAAAATTTAGTTAACAGCTAAACACCCTAAACAACCTGGCTTCAATCTAAAGCCCCGGCAGA-GTTGAAGCTGCTTCTTTGAACTTGCAATTCAACGTGAAAAATCACTTCGGAGCTTGGCAAAAAGAGGTTTCACCTCTGTCCTTAGATTTACAGTCTAATGCTTTA-CTCAGCCACTTTACCACAAAAAAGGAAGGAATCGAACCCTCTAAAACCGGTTTCAAGCCAGCCCCATAACCTTTATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATCACAGGTCCAAACCCCGTATATCTTATCACTGTAGAGCTAGACCAGCATTAACCTTTTAAGTTAAAGACTAAGAGAACTACCGCCTCTTTACAGTGA"}};



	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		
		System.setErr(new PrintStream(new OutputStream() {
		    @Override
			public void write(int b) {
		    }
		}));
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
	public void testTreeLikelihoodCustomData(){

		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/Stage0/";
		
		if( Files.exists(Paths.get(dataDir)) ){

			String trueAlignmentFile = "121101_true_seqs.fasta";
			String truePhylogenyFile = "121101_true_tree.newick";
//			String shortReadFile = "121101_short_reads.fasta";
//			String refSeqFile = "121101_ref.fasta";
			
			DataImporter dataImporter = new DataImporter(dataDir);
			SimpleAlignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
			Tree truePhylogeny = dataImporter.importTree(truePhylogenyFile);
//			Sequences shortReads = dataImporter.importSequence(shortReadFile);
//			Sequence refSeq = dataImporter.importRefSeq(refSeqFile);
			
			TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);
			LikelihoodCalculation li = new LikelihoodCalculation(treeModel, trueAlignment);
//			li.setTreeAndAlignment(treeModel, trueAlignment);
			assertEquals("Treelikelihood, compare to paup*", -2468.588, li.getTreeLikelihood(), 1e-3);
			
		}
		else{
			fail("change hardcoded dataDir: "+dataDir);
		}
	}
	

	@Test
	public void testCoalescentLikelihoodCustomData(){

		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/Stage0/";
		
		if( Files.exists(Paths.get(dataDir)) ){

			String trueAlignmentFile = "121101_true_seqs.fasta";
			String truePhylogenyFile = "121101_true_tree.newick";
			String shortReadFile = "121101_short_reads.fasta";
			String refSeqFile = "121101_ref.fasta";
			
			DataImporter dataImporter = new DataImporter(dataDir);
			SimpleAlignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
			Tree truePhylogeny = dataImporter.importTree(truePhylogenyFile);
//			Sequences shortReads = dataImporter.importSequence(shortReadFile);
//			Sequence refSeq = dataImporter.importRefSeq(refSeqFile);
			
			TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);
			LikelihoodCalculation li = new LikelihoodCalculation(treeModel, trueAlignment);
			li.setPopSize(3000, 0, 30000);

//			assertEquals("coalescent likelihood", -2468.588, li.getCoalescentLikelhood(), 1e-3);
			
		}
		else{
			fail("change hardcoded dataDir: "+dataDir);
		}
	}
	

	@Test
	public void testShortReadLikelihoodCustomData(){
		//Note: Take way too long to test this on 1000 short reads
		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/Stage0/";
		
		if( Files.exists(Paths.get(dataDir)) ){

			String trueAlignmentFile = "121101_true_seqs.fasta";
			String truePhylogenyFile = "121101_true_tree.newick";
			String shortReadFile = "121101_short_reads.fasta";
//			String refSeqFile = "121101_ref.fasta";
			
			DataImporter dataImporter = new DataImporter(dataDir);
			SimpleAlignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
			Tree truePhylogeny = dataImporter.importTree(truePhylogenyFile);
			Sequences shortReads = dataImporter.importSequence(shortReadFile);
//			Sequence refSeq = dataImporter.importRefSeq(refSeqFile);
			
			TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);
			
			shortReadFile = "121101_short_reads_1.fasta";
			shortReads = dataImporter.importSequence(shortReadFile);
			ShortReadLikelihood srl = new ShortReadLikelihood(shortReads, trueAlignment);
			LikelihoodCalculation li = new LikelihoodCalculation(treeModel, trueAlignment, shortReads);
			li.setShortReads(shortReads);
			double expected = -0.994729505252;
			assertEquals("Compare to python version, 1 read", expected, srl.getLogLikelihood(), 1e-10);
			assertEquals("Compare to python version, 1 read", expected, li.getShortReadLikelihood(), 1e-10);
			
			shortReadFile = "121101_short_reads_10.fasta";
			shortReads = dataImporter.importSequence(shortReadFile);
			srl = new ShortReadLikelihood(shortReads, trueAlignment);
			li = new LikelihoodCalculation(treeModel, trueAlignment, shortReads);
			expected = -11.855233700763;
			assertEquals("Compare to python version, 10 reads", expected, srl.getLogLikelihood(), 1e-10);
			assertEquals("Compare to python version, 10 read", expected, li.getShortReadLikelihood(), 1e-10);

			shortReadFile = "121101_short_reads_100.fasta";
			shortReads = dataImporter.importSequence(shortReadFile);
			srl = new ShortReadLikelihood(shortReads, trueAlignment);
			li = new LikelihoodCalculation(treeModel, trueAlignment, shortReads);
			expected = -106.642313524702;
			//TODO: Add this back for complete test unit, ~ +15s to run time
//			assertEquals("Compare to python version, 100 reads", expected, srl.getLogLikelihood(), 1e-10);
//			assertEquals("Compare to python version, 100 read", expected, li.getShortReadLikelihood(), 1e-10);
		
		
		}
		else{
			fail("change hardcoded dataDir: "+dataDir);
		}
	}
	
	
	@Test
    public void testNewickTree() {
//        System.out.println("\nTest Simple Node to convert Newick Tree:");
        String expectedNewickTree = "((((human:0.024003,(chimp:0.010772,bonobo:0.010772):0.013231):0.012035," +
                "gorilla:0.036038):0.033087,orangutan:0.069125):0.030457,siamang:0.099582);";
        SimpleAlignment alignment = createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
	    TreeModel treeModel = createPrimateTreeModel (alignment);
	    
	    assertEquals(expectedNewickTree, Tree.Utils.newick(treeModel, 6));
        
	   
	    
    }


	@Test
	public void testTreeLikelihoodJC69() {

	    
	    SimpleAlignment alignment = createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
	    TreeModel treeModel = createPrimateTreeModel (alignment);
	    
		LikelihoodCalculation lic = new LikelihoodCalculation();
	    lic.setTreeAndAlignment(treeModel, alignment);
	    
	    double actual = lic.getTreeLikelihood();
	    double expected = -1992.20564;
	    
	    assertEquals("treeLikelihoodJC69", expected , actual, 1e-5);
	    
	    NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);
	    format.setMaximumFractionDigits(5);
        assertEquals("treeLikelihoodJC69", format.format(expected), format.format(actual));

		    
	}


	@Test
	public void testCoalescent() throws Exception {
		/*
			#R Code
data<- c(0,0.030457, 0.033087, 0.012035, 0.013231, 0.010772)
N<- 3000

coalescentLikelihood<- function(i, N, t){
	l = (1/N)*exp( -i*(i-1)/(2*N)*t    )
	return(l)
}
logCoalescentLikelihood<- function(i, N, t){
	l = log(1/N) + ( -i*(i-1)/(2*N)*t    )
	return(l)
}

l=1
for(i in 2:length(data)){
	l = l*coalescentLikelihood(i,N,data[i] )
}

logL=0
for(i in 2:length(data)){
	logL = logL+ logCoalescentLikelihood(i,N,data[i] ) 
}


log(l)
logL


		 */		
		SimpleAlignment alignment = createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
		TreeModel treeModel = createPrimateTreeModel (alignment);
	    
    	Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000,3000,3000);
    	ConstantPopulationModel constantModel = new ConstantPopulationModel(popSize, Units.Type.DAYS);//createRandomInitialTree(popSize);

    	CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel, null, new ArrayList<TaxonList>(), constantModel);
//    	coalescent.setId("coalescent");

		LikelihoodCalculation li = new LikelihoodCalculation(treeModel, alignment, popSize);
		double expected = -40.03200311091789;
		assertEquals("coalescent", expected, li.getCoalescentLikelhood(), 1e-10);
		assertEquals("coalescent", expected, coalescent.getLogLikelihood(), 1e-10);
		

		popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 5000,5000,5000);
		li.setPopSize(5000,0,500000);
		expected =  -42.5860651207;
		assertEquals("coalescent", expected, li.getCoalescentLikelhood(), 1e-10);

    }


    @Test
    public void testStrictClock2() throws Exception {
    	SimpleAlignment alignment = createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    	TreeModel treeModel = createPrimateTreeModel (alignment);
    	Sequences reads = new Sequences();

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

    	TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null,
    			false, false, true, false, false);
    	treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);

    	//shortReadLikelihood
    	DataImporter dataImporter = new DataImporter("/home/sw167/workspace/ABI/jar/");
		String shortReadFile = "121101_short_reads_1.fasta";
		Sequences shortReads = dataImporter.importSequence(shortReadFile);
		
    	Parameter ShortRead = new Parameter.Default("SHORTREAD", 1.0, 0, 100.0);
//    	ShortReadLikelihood shortReadLikelihood = new ShortReadLikelihood(shortReads, alignment);
    	ShortReadLikelihood shortReadLikelihood = new ShortReadLikelihood(shortReads, alignment, ShortRead);
    	
    	// Operators
    	OperatorSchedule schedule = new SimpleOperatorSchedule();

    	MCMCOperator operator = new ScaleOperator(kappa, 0.75);
    	schedule.addOperator(operator);

    	operator = new ScaleOperator(rateParameter, 0.75);
    	operator.setWeight(3.0);
    	schedule.addOperator(operator);

    	Parameter allInternalHeights = treeModel.createNodeHeightsParameter(true, true, false);
    	operator = new UpDownOperator(new Scalable[]{new Scalable.Default(rateParameter)},
    			new Scalable[] {new Scalable.Default(allInternalHeights)}, 0.75, 3.0, CoercionMode.COERCION_ON);
//    	schedule.addOperator(operator);

    	operator = new ScaleOperator(popSize, 0.75);
    	operator.setWeight(3.0);
    	schedule.addOperator(operator);

    	Parameter rootHeight = treeModel.getRootHeightParameter();
    	String TREE_HEIGHT = TreeModel.TREE_MODEL + "." + TreeModelParser.ROOT_HEIGHT;
    	rootHeight.setId(TREE_HEIGHT);
    	operator = new ScaleOperator(rootHeight, 0.75);
    	operator.setWeight(3.0);
//    	schedule.addOperator(operator);

    	Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);
    	operator = new UniformOperator(internalHeights, 30.0);
//    	schedule.addOperator(operator);

    	operator = new SubtreeSlideOperator(treeModel, 15.0, 1.0, true, false, false, false, CoercionMode.COERCION_ON);
//    	schedule.addOperator(operator);

    	operator = new ExchangeOperator(ExchangeOperator.NARROW, treeModel, 15.0);
    	//         operator.doOperation();
//    	schedule.addOperator(operator);

    	operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 3.0);
    	//         operator.doOperation();
//    	schedule.addOperator(operator);

    	operator = new WilsonBalding(treeModel, 3.0);
    	//         operator.doOperation();
//    	schedule.addOperator(operator);

    	
    	
    	//test new operator
    	operator = new AlignmentOperator(ShortRead, alignment, 1);
    	operator.setWeight(10);
    	schedule.addOperator(operator);

    	
    	
    	//CompoundLikelihood
    	List<Likelihood> likelihoods = new ArrayList<Likelihood>();        
    	likelihoods.add(coalescent);
    	Likelihood prior = new CompoundLikelihood(0, likelihoods);
    	prior.setId(CompoundLikelihoodParser.PRIOR);

    	likelihoods.clear();
    	likelihoods.add(treeLikelihood);
    	Likelihood likelihood = new CompoundLikelihood(-1, likelihoods);

    	likelihoods.clear();
    	likelihoods.add(shortReadLikelihood);
    	Likelihood shortReadLikelihoodCompound = new CompoundLikelihood(-1, likelihoods);
    	shortReadLikelihoodCompound.setId("ShortReadLikelihood");
    	
    	likelihoods.clear();
    	likelihoods.add(prior);
    	likelihoods.add(likelihood);
    	likelihoods.add(shortReadLikelihood);
    	Likelihood posterior = new CompoundLikelihood(0, likelihoods);
    	posterior.setId(CompoundLikelihoodParser.POSTERIOR);

    	//         System.out.println(posterior.getId());
    	//         System.out.println(posterior.getLogLikelihood());
    	//         System.out.println(coalescent.getLogLikelihood());
    	//         System.out.println(likelihood.getLogLikelihood());
    	//         System.out.println(posterior.getModel().getModelName());
    	//         System.out.println(posterior.prettyName());
    	//         
    	//         CompoundLikelihood l2 = new CompoundLikelihood(0, likelihoods);
    	//         System.out.println(l2.getDiagnosis());
    	//         System.out.println(l2.getLikelihoodCount());
    	//         System.out.println(l2.getReport());
    	//         System.out.println(l2.toString());
    	//         System.out.println();

    	ArrayLogFormatter formatter = new ArrayLogFormatter(false);
    	
    	int lengthScaler = 1;
    	MCLogger[] loggers = new MCLogger[2];
    	loggers[0] = new MCLogger(formatter, lengthScaler*1, false);
    	loggers[0].add(posterior);
    	loggers[0].add(treeLikelihood);
    	loggers[0].add(rootHeight);
    	loggers[0].add(rateParameter);
    	loggers[0].add(popSize);
    	loggers[0].add(kappa);
    	loggers[0].add(coalescent);

    	loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), lengthScaler*1, false);
    	loggers[1].add(posterior);
    	loggers[1].add(treeLikelihood);
//    	loggers[1].add(coalescent);
    	loggers[1].add(prior);
    	loggers[1].add(shortReadLikelihoodCompound);
//    	loggers[1].add(rootHeight);
    	loggers[1].add(rateParameter);
//    	loggers[1].add(mu);
    	loggers[1].add(popSize);

    	//////////////////
//    	fail("fail here to stop MCMC runs");

    	// MCMC
    	MCMC mcmc = new MCMC("mcmc1");
    	MCMCOptions options = new MCMCOptions();
    	options.setChainLength(1000000);
    	options.setChainLength(lengthScaler*10);
//    	options.setUseCoercion(true); // autoOptimize = true
//    	options.setCoercionDelay(lengthScaler*5);
//    	options.setTemperature(1.0);
//    	options.setFullEvaluationCount(lengthScaler*2);

    	mcmc.setShowOperatorAnalysis(true);
    	mcmc.init(options, posterior, schedule, loggers);
    	mcmc.run();

    	// time
    	System.out.println(mcmc.getTimer().toString());

    	// Tracer
    	List<Trace> traces = formatter.getTraces();
    	ArrayTraceList traceList = new ArrayTraceList("RandomLocalClockTest", traces, 0);

    	for (int i = 1; i < traces.size(); i++) {
    		traceList.analyseTrace(i);
    	}


//        <expectation name="posterior" value="-3928.71"/>
//        <expectation name="clock.rate" value="8.04835E-4"/>
//        <expectation name="constant.popSize" value="37.3762"/>
//        <expectation name="hky.kappa" value="18.2782"/>
//        <expectation name="treeModel.rootHeight" value="69.0580"/>
//        <expectation name="treeLikelihood" value="-3856.59"/>
//        <expectation name="coalescent" value="-72.1285"/>

    	TraceCorrelation likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(CompoundLikelihoodParser.POSTERIOR));
//        assertExpectation(CompoundLikelihoodParser.POSTERIOR, likelihoodStats, -3928.71);
//
//        likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TreeLikelihoodParser.TREE_LIKELIHOOD));
//        assertExpectation(TreeLikelihoodParser.TREE_LIKELIHOOD, likelihoodStats, -3856.59);
//
//        TraceCorrelation treeHeightStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TREE_HEIGHT));
//        assertExpectation(TREE_HEIGHT, treeHeightStats, 69.0580);
//
//        TraceCorrelation kappaStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(HKYParser.KAPPA));
//        assertExpectation(HKYParser.KAPPA, kappaStats, 18.2782);
//
//        TraceCorrelation rateStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(StrictClockBranchRates.RATE));
//        assertExpectation(StrictClockBranchRates.RATE, rateStats, 8.04835E-4);        
//
//        TraceCorrelation popStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(ConstantPopulationModelParser.POPULATION_SIZE));
//        assertExpectation(ConstantPopulationModelParser.POPULATION_SIZE, popStats, 37.3762);
//
//        TraceCorrelation coalescentStats = traceList.getCorrelationStatistics(traceList.getTraceIndex("coalescent"));
//        assertExpectation("coalescent", coalescentStats, -72.1285);
    }

    //************************** data ****************************
	
    private OperatorSchedule setupDefaultSchedule(TreeModel treeModel, Parameter kappa, Parameter rateParameter, Parameter popSize, Parameter rootHeight) {

        OperatorSchedule schedule = new SimpleOperatorSchedule();
        
        MCMCOperator operator = new ScaleOperator(kappa, 0.75);
        operator.setWeight(1.0);

		schedule.addOperator(operator);

        operator = new ScaleOperator(rateParameter, 0.75);
        operator.setWeight(3.0);
        schedule.addOperator(operator);

        Parameter allInternalHeights = treeModel.createNodeHeightsParameter(true, true, false);
        operator = new UpDownOperator(new Scalable[]{new Scalable.Default(rateParameter)},
                new Scalable[] {new Scalable.Default(allInternalHeights)}, 0.75, 3.0, CoercionMode.COERCION_ON);
        schedule.addOperator(operator);

        operator = new ScaleOperator(popSize, 0.75);
        operator.setWeight(3.0);
        schedule.addOperator(operator);

        
        operator = new ScaleOperator(rootHeight, 0.75);
        operator.setWeight(3.0);
        schedule.addOperator(operator);

        Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);
        operator = new UniformOperator(internalHeights, 30.0);
        schedule.addOperator(operator);

        operator = new SubtreeSlideOperator(treeModel, 15.0, 1.0, true, false, false, false, CoercionMode.COERCION_ON);
        schedule.addOperator(operator);

        operator = new ExchangeOperator(ExchangeOperator.NARROW, treeModel, 15.0);
//        operator.doOperation();
        schedule.addOperator(operator);

        operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 3.0);
//        operator.doOperation();
        schedule.addOperator(operator);

        operator = new WilsonBalding(treeModel, 3.0);
//        operator.doOperation();
        schedule.addOperator(operator);
		return schedule;
	}

    private static SimpleAlignment createAlignment(Object[][] taxa_sequence, DataType dataType) {
	
	        SimpleAlignment alignment = new SimpleAlignment();
	        alignment.setDataType(dataType);
	//        alignment.setDataType(Nucleotides.INSTANCE);
	        
	        Taxon[] taxa = new Taxon[taxa_sequence[0].length]; // 6, 17
	//        System.out.println("Taxon len = " + taxa_sequence[0].length);
	//        System.out.println("Alignment len = " + taxa_sequence[1].length);
	        if (taxa_sequence.length > 2) System.out.println("Date len = " + taxa_sequence[2].length);                          
	
	        for (int i=0; i < taxa_sequence[0].length; i++) {
	            taxa[i] = new Taxon(taxa_sequence[0][i].toString());
	
	            if (taxa_sequence.length > 2) {
	                Date date = new Date((Double) taxa_sequence[2][i], Units.Type.YEARS, (Boolean) taxa_sequence[3][0]);
	                taxa[i].setDate(date);
	            }
	
	            //taxonList.addTaxon(taxon);
	            Sequence sequence = new Sequence(taxa_sequence[1][i].toString());
	            sequence.setTaxon(taxa[i]);
	            sequence.setDataType(dataType);
	
	            alignment.addSequence(sequence);
	        }
	        return alignment;
	    }

	private static TreeModel createPrimateTreeModel (SimpleAlignment alignment) {
	    	
	    	Taxon[] taxa = alignment.asList().toArray(new Taxon[0]);
	    	
	        SimpleNode[] nodes = new SimpleNode[10];
	        for (int n=0; n < 10; n++) {
	            nodes[n] = new SimpleNode();
	        }
	
	//        nodes[0].setHeight(0);
	        nodes[0].setTaxon(taxa[0]); // human
	
	        nodes[1].setTaxon(taxa[1]); // chimp
	
	        nodes[2].setTaxon(taxa[2]); // bonobo
	
	        nodes[3].setHeight(0.010772);
	        nodes[3].addChild(nodes[1]);
	        nodes[3].addChild(nodes[2]);
	
	        nodes[4].setHeight(0.024003);
	        nodes[4].addChild(nodes[0]);
	        nodes[4].addChild(nodes[3]);
	
	        nodes[5].setTaxon(taxa[3]); // gorilla
	
	        nodes[6].setHeight(0.036038);
	        nodes[6].addChild(nodes[4]);
	        nodes[6].addChild(nodes[5]);
	
	        nodes[7].setTaxon(taxa[4]); // orangutan
	
	        nodes[8].setHeight(0.069125);
	        nodes[8].addChild(nodes[6]);
	        nodes[8].addChild(nodes[7]);
	
	        nodes[9].setTaxon(taxa[5]); // siamang
	
	        SimpleNode root = new SimpleNode();
	        root.setHeight(0.099582);
	        root.addChild(nodes[8]);
	        root.addChild(nodes[9]);
	
	        Tree tree = new SimpleTree(root);
	        tree.setUnits(Units.Type.YEARS);
	        
	        return new TreeModel(tree); //treeModel
	    }

	protected ConstantPopulationModel createRandomInitialTree(Parameter popSize) {        
        ConstantPopulationModel startingTree = new ConstantPopulationModel(popSize, Units.Type.DAYS);
        ConstantPopulation constant = (ConstantPopulation) startingTree.getDemographicFunction();

//        createTreeModel(constant);

        return startingTree;
    }

//    protected void createTreeModel (ConstantPopulation constant) {
//        CoalescentSimulator simulator = new CoalescentSimulator();
//        Tree tree = simulator.simulateTree(alignment, constant);
//        treeModel = new TreeModel(tree);//treeModel
//    }
	
    protected void assertExpectation(String name, TraceCorrelation stats, double v) {
        double mean = stats.getMean();
        double stderr = stats.getStdErrorOfMean();
        double upper = mean + 2 * stderr;
        double lower = mean - 2 * stderr;

        assertTrue("Expected " + name + " is " + v + " but got " + mean + " +/- " + stderr,
                upper > v && lower < v);
    }

    
}
