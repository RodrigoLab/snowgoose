package test;

import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.Tree;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
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
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.operators.UniformOperator;
import dr.inference.trace.ArrayTraceList;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceCorrelation;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.math.MathUtils;
/**
 * @author Walter Xie
 * convert testStrictClock.xml in the folder /example
 */
public class CopyStrictClockTest  {


    protected static final String TREE_HEIGHT = TreeModel.TREE_MODEL + "." + TreeModelParser.ROOT_HEIGHT;

    protected TreeModel treeModel;
    protected SimpleAlignment alignment;
    protected Taxon[] taxa;

//    public CopyStrictClockTest(String name) {
////        super(name);
//    }


	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

    @Before
    public void setUp() throws Exception {
//        super.setUp();

        MathUtils.setSeed(666);

        createAlignment(DENGUE4_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }
    
    

    @Test
    public void testStrictClock() throws Exception {
        Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 380.0, 0, 38000.0);
        ConstantPopulationModel constantModel = createRandomInitialTree(popSize);

        CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel, null, new ArrayList<TaxonList>(), constantModel);
        coalescent.setId("coalescent");

        // clock model
        Parameter rateParameter =  new Parameter.Default(StrictClockBranchRates.RATE, 2.3E-5, 0, 100.0);

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

        // Operators
        OperatorSchedule schedule = new SimpleOperatorSchedule();

        MCMCOperator operator = new ScaleOperator(kappa, 0.75);
        operator.setWeight(1.0);
        schedule.addOperator(operator);

//        operator = new ScaleOperator(rateParameter, 0.75);
//        operator.setWeight(3.0);
//        schedule.addOperator(operator);

//        Parameter allInternalHeights = treeModel.createNodeHeightsParameter(true, true, false);
//        operator = new UpDownOperator(new Scalable[]{new Scalable.Default(rateParameter)},
//                new Scalable[] {new Scalable.Default(allInternalHeights)}, 0.75, 3.0, CoercionMode.COERCION_ON);
//        schedule.addOperator(operator);

//        operator = new ScaleOperator(popSize, 0.75);
//        operator.setWeight(3.0);
//        schedule.addOperator(operator);

//        Parameter rootHeight = treeModel.getRootHeightParameter();
//        rootHeight.setId(TREE_HEIGHT);
//        operator = new ScaleOperator(rootHeight, 0.75);
//        operator.setWeight(3.0);
//        schedule.addOperator(operator);

        Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);
        operator = new UniformOperator(internalHeights, 30.0);
        schedule.addOperator(operator);

        operator = new SubtreeSlideOperator(treeModel, 15.0, 1.0, true, false, false, false, CoercionMode.COERCION_ON);
        schedule.addOperator(operator);

//        operator = new ExchangeOperator(ExchangeOperator.NARROW, treeModel, 15.0);
////        operator.doOperation();
//        schedule.addOperator(operator);
//
//        operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 3.0);
////        operator.doOperation();
//        schedule.addOperator(operator);

        operator = new WilsonBalding(treeModel, 3.0);
//        operator.doOperation();
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
        likelihoods.add(prior);
        likelihoods.add(likelihood);
        Likelihood posterior = new CompoundLikelihood(0, likelihoods);
        posterior.setId(CompoundLikelihoodParser.POSTERIOR);

        // Log
        ArrayLogFormatter formatter = new ArrayLogFormatter(false);

        MCLogger[] loggers = new MCLogger[2];
        loggers[0] = new MCLogger(formatter, 500, false);
        loggers[0].add(posterior);
        loggers[0].add(treeLikelihood);
//        loggers[0].add(rootHeight);
        loggers[0].add(rateParameter);
        loggers[0].add(popSize);
        loggers[0].add(kappa);
        loggers[0].add(coalescent);

        loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), 10000, false);
        loggers[1].add(posterior);
        loggers[1].add(treeLikelihood);
//        loggers[1].add(rootHeight);
        loggers[1].add(rateParameter);
        loggers[1].add(coalescent);

        // MCMC
        MCMC mcmc = new MCMC("mcmc1");
        MCMCOptions options = new MCMCOptions(1000000);

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
        assertExpectation(CompoundLikelihoodParser.POSTERIOR, likelihoodStats, -3928.71);

        likelihoodStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TreeLikelihoodParser.TREE_LIKELIHOOD));
        assertExpectation(TreeLikelihoodParser.TREE_LIKELIHOOD, likelihoodStats, -3856.59);

        TraceCorrelation treeHeightStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TREE_HEIGHT));
        assertExpectation(TREE_HEIGHT, treeHeightStats, 69.0580);

        TraceCorrelation kappaStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(HKYParser.KAPPA));
        assertExpectation(HKYParser.KAPPA, kappaStats, 18.2782);

        TraceCorrelation rateStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(StrictClockBranchRates.RATE));
        assertExpectation(StrictClockBranchRates.RATE, rateStats, 8.04835E-4);        

        TraceCorrelation popStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(ConstantPopulationModelParser.POPULATION_SIZE));
        assertExpectation(ConstantPopulationModelParser.POPULATION_SIZE, popStats, 37.3762);

        TraceCorrelation coalescentStats = traceList.getCorrelationStatistics(traceList.getTraceIndex("coalescent"));
        assertExpectation("coalescent", coalescentStats, -72.1285);
    }


    protected void createAlignment(Object[][] taxa_sequence, DataType dataType) {

        alignment = new SimpleAlignment();
        alignment.setDataType(dataType);
//        alignment.setDataType(Nucleotides.INSTANCE);

        taxa = new Taxon[taxa_sequence[0].length]; // 6, 17
        System.out.println("Taxon len = " + taxa_sequence[0].length);
        System.out.println("Alignment len = " + taxa_sequence[1].length);
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
    }
    protected void assertExpectation(String name, TraceCorrelation stats, double v) {
        double mean = stats.getMean();
        double stderr = stats.getStdErrorOfMean();
        double upper = mean + 2 * stderr;
        double lower = mean - 2 * stderr;

        assertTrue("Expected " + name + " is " + v + " but got " + mean + " +/- " + stderr,
                upper > v && lower < v);
    }

    protected static final Object[][] DENGUE4_TAXON_SEQUENCE = {{"D4Brazi82", "D4ElSal83", "D4ElSal94", "D4Indon76",
            "D4Indon77", "D4Mexico84", "D4NewCal81", "D4Philip64", "D4Philip56", "D4Philip84", "D4PRico86", "D4SLanka78",
            "D4Tahiti79", "D4Tahiti85", "D4Thai63", "D4Thai78", "D4Thai84"}, // 17
          {"ATGCGATGCGTAGGAGTAGGAAACAGAGACTTTGTGGAAGGAGTCTCAGGTGGAGCATGGGTCGACCTGGTGCTAGAACATGGAGGATGCGTCACAACCATGGCCCAGGGAAAACCAACCTTGGATTTTGAACTGACCAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCAATATCAAACATAACTACGGCAACAAGATGTCCAACGCAAGGAGAGCCTTATCTGAAAGAGGAACAGGACCAACAGTACATTTGCCGGAGAGATGTGGTAGACAGAGGGTGGGGCAATGGCTGTGGCTTGTTTGGAAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGTTCGGGGAAGATAACAGGCAATTTGGTCCAAATTGAGAACCTTGAATACACAGTGGTTGTAACAGTCCACAATGGAGACACCCATGCAGTAGGAAATGACACATCCAATCATGGAGTTACAGCCATGATAACTCCCAGGTCACCATCGGTGGAAGTCAAATTGCCGGACTATGGAGAACTAACACTCGATTGTGAACCCAGGTCTGGAATTGACTTTAATGAGATGATTCTGATGAAAATGAAAAAGAAAACATGGCTCGTGCATAAGCAATGGTTTTTGGATCTGCCTCTTCCATGGACAGCAGGAGCAGACACATCAGAGGTTCACTGGAATTACAAAGAGAGAATGGTGACATTTAAGGTTCCTCATGCCAAGAGACAGGATGTGACAGTGCTGGGATCTCAGGAAGGAGCCATGCATTCTGCCCTCGCTGGAGCCACAGAAGTGGACTCCGGTGATGGAAATCACATGTTTGCAGGACATCTCAAGTGCAAAGTCCGTATGGAGAAATTGAGAATCAAGGGAATGTCATACACGATGTGTTCAGGAAAGTTTTCAATTGACAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTGAAAGTCAAGTATGAAGGTGCTGGAGCTCCGTGTAAAGTCCCCATAGAGATAAGAGATGTAAACAAGGAAAAAGTGGTTGGGCGTATCATCTCATCCACCCCTTTGGCTGAGAATACCAACAGTGTAACCAACATAGAATTAGAACCCCCCTTTGGGGACAGCTACATAGTGATAGGTGTTGGAAACAGCGCATTAACACTCCATTGGTTCAGGAAAGGGAGTTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGTGCAAAACGAATGGCCATTCTAGGTGAAACAGCTTGGGATTTTGGTTCCGTTGGTGGATTGTTCACATCATTGGGAAAGGCTGTGCACCAGGTTTTTGGAAGTGTGTATACAACCATGTTTGGAGGAGTCTCATGGATGATTAGAATCCTAATTGGGTTCTTAGTGTTGTGGATTGGCACGAACTCAAGGAACACTTCAATGGCTATGACGTGCATAGCTGTTGGAGGAATCACTCTGTTTCTGGGCTTCACAGTTCAAGCA",
           "ATGCGATGCGTAGGAGTAGGAAACAGAGACTTTGTGGAAGGAGTCTCAGGTGGAGCATGGGTCGACCTGGTGCTAGAACATGGAGGATGCGTCACAACCATGGCCCAGGGAAAACCAACCTTGGATTTTGAACTGACTAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCAATATCAAACATAACTACGGCAACAAGATGTCCAACGCAAGGAGAGCCTTATCTGAAAGAGGAACAGGACCAACAGTACATTTGCCGGAGAGATGTGGTAGACAGAGGGTGGGGCAATGGCTGTGGCTTGTTTGGAAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGTTCGGGGAAGATAACAGGCAATTTGGTCCAAATTGAGAACCTTGAATACACAGTGGTTGTAACAGTCCACAATGGAGACACCCATGCAGTAGGAAATGACACATCCAATCATGGAGTTACAGCCATGATAACTCCCAGGTCACCATCGGTGGAAGTCAAATTGCCGGACTATGGAGAACTAACACTCGATTGTGAACCCAGGTCTGGAATTGACTTTAATGAGATGATTCTGATGAAAATGAAAAAGAAAACATGGCTCGTGCATAAGCAATGGTTTTTGGATCTGCCTCTTCCATGGACAGCAGGAGCAGACACATCAGAGGTTCACTGGAATTACAAAGAGAGAATGGTGACATTTAAGGTTCCTCATGCCAAGAGACAGGATGTGACAGTGCTGGGATCTCAGGAAGGAGCCATGCATTCTGCCCTCGCTGGAGCCACAGAAGTGGACTCCGGTGATGGAAATCATATGTTTGCAGGACATCTCAAGTGCAAAGTCCGTATGGAGAAATTGAGAATCAAGGGAATGTCATACACGATGTGTTCAGGAAAGTTTTCAATTGACAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTGAAAGTCAAGTATGAAGGTGCTGGAGCTCCGTGTAAAGTCCCCATAGAGATAAGAGATGTAAACAAGGAAAAAGTGGTTGGGCGTATCATCTCATCCACCCCTTTGGCTGAGAATACCAACAGTGTAACCAACATAGAATTAGAACCCCCCTTTGGGGACAGCTACATAGTGATAGGTGTTGGAAACAGCGCATTAACACTCCATTGGTTCAGGAAAGGGAGTTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGTGCAAAACGAATGGCCATTCTAGGTGAAACAGCTTGGGATTTTGGTTCCGTTGGTGGACTGTTCACATCATTGGGAAAGGCTGTGCACCAGGTTTTTGGAAGTGTGTATACAACCATGTTTGGAGGAGTCTCATGGATGATTAGAATCCTAATTGGGTTCTTAGTGTTGTGGATTGGCACGAACTCAAGGAACACTTCAATGGCTATGACGTGCATAGCTGTTGGAGGAATCACTCTGTTTCTGGGCTTCACAGTTCAAGCA",
           "ATGCGATGCGTAGGAGTAGGAAACAGAGACTTTGTGGAAGGAGTCTCAGGTGGAGCATGGGTCGACCTGGTGCTAGAACATGGAGGATGCGTCACAACCATAGCCCAGGGAAAACCAACCTTGGATTTTGAATTGACTAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCAATATCAAACATAACTACGGCAACAAGATGTCCAACGCAAGGAGAGCCTTATCTGAAAGAGGAACAGGACCAACAGTACATTTGCCGGAGAGATGTGGTAGACAGAGGGTGGGGGAATGGCTGTGGCTTGCTTGGAAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGTTCGGGGAAGATAACAGGCAATTTGGTCCAAATTGAGAACCTTGAATACACAGTGGTTGTAACAGTCCACAATGGAGATACCCATGCAGTAGGAAATGACACATCCAATCATGGAGTTACAGCCACGATAACTCCCAGGTCACCATCGGTGGAAGTCAAATTGCCGGACTATGGAGAACTAACACTCGATTGTGAACCCAGATCTGGAATTGATTTTAATGAGATGATTCTGATGAAAATGAAAAAGAAAACATGGCTCGTGCATAAGCAATGGTTTTTGGATCTGCCTCTTCCATGGACAGCAGGAGCAGACACATCAGAGGTTCACTGGAATTACAAAGAGAGAATGGTGACATTCAAGGTTCCTCATGCCAAGAGACAGGATGTGACAGTGCTGGGATCTCAGGAAGGAGCCATGCATTCTGCCCTCGCTGGAGCCACAGAAGTGGACTCCGGTGATGGAAATCACATGTTTGCAGGACATCTCAAGTGCAAAGTCCGCATGGAGAAATTGAGAATCAAGGGAATGTCATACACGATGTGTTCAGGAAAGTTTTCAATTGATAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTGAAAGTCAAGTATGAAGGTGCTGGAGCTCCGTGTAAAGTCCCCATAGAGATAAGAGATGTAAACAAGGAAAAAGTGGTTGGGCGTATCATCTCATCCACCCCTTTGGCTGAGAATACCAACAGTGTAACCAACATAGAATTAGAACCCCCCTTTGGGGACAGCTACATAGTGATAGGTGTCGGAAACAGCGCATTAACACTCCATTGGTTCAGGAAAGGGAGTTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGTGCAAAACGAATGGCCATTCTAGGTGAAACAGCTTGGGATTTTGGTTCCGTTGGTGGACTGTTCACATCATTGGGAAAGGCTGTGCACCAGGTTTTTGGAAGTGTGTACACAACCATGTTTGGAGGAGTCTCATGGATGATTAGAATCCTAATTGGGTTCTTAGTGTTATGGATTGGCACGAACTCAAGGAACACTTCAATGGCTATGACGTGCATAGCTGTTGGAGGAATCACTCTGTTTCTGGGCTTCACAGCTCAAGCA",
           "ATGCGATGCGTAGGAGTAGGAAACAGAGACTTTGTGGAAGGAGTCTCAGGTGGAGCATGGGTCGATCTGGTGCTAGAACATGGAGGATGCGTCACAACCATGGCCCAGGGAAAACCAACCTTGGATTTTGAACTGACTAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCAATATCAAACATAACCACGGCAACAAGATGTCCAACGCAAGGAGAGCCTTATCTAAAAGAGGAACAAGACCAACAGTACATTTGCCGGAGAGATGTGGTAGACAGAGGGTGGGGCAATGGCTGTGGCTTGTTTGGAAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGTTCGGGGAAGATAACAGGCAATTTGGTCCAAATTGAGAACCTTGAATACACAGTGGTTGTAACAGTCCACAATGGAGACACCCATGCAGTAGGAAATGACACATCCAATCATGGAGTTACAGCCACGATAACTCCCAGGTCACCATCGGTGGAAGTCAAATTGCCGGACTATGGAGAACTAACACTCGATTGTGAACCCAGGTCTGGAATTGACTTTAATGAGATGATTCTGATGAAAATGAAAAAGAAAACATGGCTTGTGCATAAGCAATGGTTTTTGGATCTACCTCTACCATGGACAGCAGGAGCAGACACATCAGAGGTTCACTGGAATTACAAAGAGAGAATGGTGACATTTAAGGTTCCTCATGCCAAGAGACAGGATGTGACAGTGCTGGGATCTCAGGAAGGAGCCATGCATTCTGCCCTCGCTGGAGCCACAGAAGTGGACTCCGGTGATGGAAATCACATGTTTGCAGGACATCTCAAGTGCAAAGTCCGTATGGAGAAATTGAGAATCAAGGGAATGTCATACACGATGTGTCCAGGAAAGTTCTCAATTGACAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTGAAAGTCAAGTATGAAGGTGCTGGAGCTCCGTGTAAAGTCCCCATAGAGATAAGAGATGTGAACAAGGAAAAAGTGGTTGGGCGTATCATCTCATCCACCCCTTTGGCTGAGAATACCAACAGTGCAACCAACATAGAGTTAGAACCCCCCTTTGGGGACAGCTACATAGTGATAGGTGTTGGAAACAGTGCATTAACACTCCATTGGTTCAGGAAAGGGAGTTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGTGCAAAACGAATGGCCATTCTAGGTGAAACAGCTTGGGATTTTGGTTCCGTTGGTGGACTGCTCACATCATTGGGAAAGGCTGTGCACCAGGTTTTTGGAAGTGTGTATACAACCATGTTTGGAGGAGTCTCATGGATGATTAGAATCCTAATTGGGTTCCTAGTGTTGTGGATTGGCACGAACTCAAGGAACACTTCAATGGCTATGACGTGCATAGCTGTTGGAGGAATCACTCTGTTTCTGGGCTTCACAGTTCAAGCA",
           "ATGCGATGCGTAGGAGTAGGAAACAGAGACTTTGTGGAAGGAGTCTCAGGTGGAGCATGGGTCGATCTGGTGCTAGAACATGGAGGATGCGTCACAACCATGGCCCAGGGAAAACCAACCTTGGATTTTGAACTGACTAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCAATATCAAACATAACCACGGCAACAAGATGTCCAACGCAAGGAGAGCCTTATCTAAAAGAGGAACAAGACCAACAGTACATTTGCCGGAGAGATGTGGTAGACAGAGGGTGGGGCAATGGCTGTGGCTTGTTTGGAAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGTTCGGGGAAGATAACAGGCAATTTGGTCCAAATTGAGAACCTTGAATACACAGTAGTTGTAACAGTCCACAATGGAGACACCCATGCAGTAGGAAATGACACATCCAACCATGGAGTTACAGCCACGATAACTCCCAGGTCACCATCGGTGGAAGTCAAATTGCCGGACTATGGAGAACTAACACTCGATTGTGAACCCAGGTCTGGAATTGACTTTAATGAGATGATTCTGATGAAAATGAAAAAGAAAACATGGCTTGTGCATAAGCAATGGTTTTTGGATCTACCTCTACCATGGACAGCAGGAGCAGACACATCAGAGGTTCACTGGAATTACAAAGAGAGAATGGTGACATTTAAGGTTCCTCATGCCAAGAGACAGGATGTGACAGTGCTGGGATCTCAGGAAGGAGCCATGCATTCTGCCCTCGCTGGAGCCACAGAAGTGGACTCCGGTGATGGAAATCACATGTTTGCAGGACATCTCAAGTGCAAAGTCCGTATGGAGAAATTGAGAATCAAGGGAATGTCATACACGATGTGTTCAGGAAAGTTCTCAATTGACAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTGAAAGTCAAGTATGAAGGTGCTGGAGCTCCGTGCAAAGTCCCCATAGAGATAAGAGATGTAAACAAGGAAAAAGTGGTTGGGCGTATCATCTCATCCACCCCTTTGGCTGAGAATACCAACAGTGTAACCAACATAGAATTAGAACCCCCCTTTGGGGACAGCTACATAGTGATAGGTGTTGGAAACAGTGCATTAACACTCCATTGGTTCAGGAAAGGGAGTTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGTGCAAAACGAATGGCCATTCTAGGTGAAACAGCTTGGGATTTTGGTTCCGTTGGTGGACTGTTCACATCATTGGGAAAGGCTGTGCACCAGGTTTTTGGAAGTGTGTATACAACCATGTTTGGAGGAGTCTCATGGATGATTAGAATCCTAATTGGCTTCTTAGTGTTGTGGATTGGCACGAACTCAAGGAACACTTCAATGGCTATGACGTGCATAGCTGTTGGAGGAATCACTCTGTTTCTGGGCTTCACAGTTCAAGCA",
           "ATGCGATGCGTAGGAGTAGGAAACAGAGACTTTGTGGAAGGAGTCTCAGGTGGAGCATGGGTCGACCTAGTGCTAGAACATGGAGGATGCGTCACAACCATGGCCCAGGGAAAACCAACCTTGGATTTTGAACTGACTAAGACAACAGCCAAGGAAGTGGCTCTGCTAAGAACCTATTGCATTGAAGCCTCAATATCAAACATAACTACGGCAACAAGATGTCCAACGCAAGGAGAGCCTTATCTGAAAGAGGAACAGGACCAACAGTACATTTGCCGGAGAGATGTGGTAGACAGAGGGTGGGGCAATGGCTGTGGCTTGTTTGGAAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGTTCGGGGAAGATAACAGGCAATTTGGTCCAAATTGAGAACCTTGAATACACAGTGGTTGTAACAGTCCACAATGGAGACACCCATGCAGTAGGAAATGACACATCCAATCATGGAGTTACAGCCATGATAACTCCCAGGTCACCATCGGTGGAAGTCAAATTGCCGGACTATGGAGAACTAACACTCGATTGTGAACCCAGGTCTGGAATTGACTTTAATGAGATGATTCTGATGAAAATGAAAAAGAAAACATGGCTCGTGCATAAGCAATGGTTTTTGGATCTGCCTCTTCCATGGACAGCAGGAGCAGACACATCAGAGGTTCACTGGAATTACAAAGAGAGAATGGTGACATTTAAGGTTCCCCATGCCAAGAGACAGGATGTGACAGTGCTGGGATCTCAGGAAGGAGCCATGCATTCTGCCCTCGCTGGAGCCACAGAAGTGGACTCCGGTGATGGAAATCACATGTTTGCAGGACATCTCAAGTGCAAAGTCCGTATGGAGAAATTGAGAATCAAGGGAATGTCATACACGATGTGTTCAGGAAAGTTTTCAATTGACAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTGAAAGTCAAGTGTGAAGGTGCTGGAGCTCCCTGTAAAGTCCCCATAGAGATAAGAGATGTAAACAAGGAAAAAGTGGTTGGGCGTATCATCTCATCCACCCCTTTGGCTGAGAATACCAACAGTGTAACCAACATAGAATTAGAACCTCCCTTTGGGGACAGCTACATAGTGATAGGTGTTGGAAACAGCGCATTAACACTCCATTGGTTCAGGAAAGGGAGTTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGTGCAAAACGAATGGCCATTCTAGGTGAAACAGCTTGGGATTTTGGTTCCGTTGGTGGACTGTTCACATCATTGGGAAAGGCTGTGCACCAGGTTTTTGGAAGTGTGTATACAACCATGTTTGGAGGAGTCTCATGGATGATTAGAATCCTAATTGGGTTCTTAGTGTTGTGGATTGGCACGAACTCAAGGAACACTTCAATGGCTATGACGTGCATAGCTGTTGGAGGAATCACTCTGTTTCTGGGCTTCACAGTTCAAGCA",
           "ATGCGATGCGTAGGAGTAGGAAACAGAGACTTTGTGGAAGGAGTCTCAGGTGGAGCATGGGTCGACCTGGTGCTAGAACATGGAGGATGCGTCACAACCATGGCCCAGGGAAAACCAACCTTGGATTTTGAACTGACTAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCAATATCAAACATAACTACGGCAACAAGATGTCCAACGCAAGGAGAGCCTTATCTGAAAGAGGAACAGGACCAACAGTACATTTGCCGGAGAGATGTGGTAGACAGAGGGTGGGGCAATGGCTGTGGCTTGTTTGGAAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGTTCGGGGAAGATAACAGGCAATTTGGTCCAAATTGAGAACCTTGAATACACAGTGGTTGTAACAGTCCACAATGGAGACACCCATGCAGTAGGAAATGACACATCCAATCATGGAGTTACAGCCATGATAACTCCCAGGTCACCATCGGTGGAAGTCAAATTGCCGGACTATGGAGAACTAACACTCGATTGTGAACCCAGGTCTGGAATTGACTTTAATGAGATGATTCTGATGAAAATGAAAAAGAAAACATGGCTCGTGCATAAGCAATGGTTTTTGGATCTGCCTCTTCCATGGACAGCAGGAGCAGACACATCAGAGGTTCACTGGAATTACAAAGAGAGAATGGTGACATTTAAGGTTCCTCATGCCAAGAGACAGGATGTGACAGTGCTGGGATCTCAGGAAGGAGCCATGCATTCTGCCCTCGCTGGAGCCACAGAAGTGGACTCCGGTGATGGAAATCACATGTTTGCAGGACATCTCAAGTGCAAAGTCCGTATGGAGAAATTGAGAATCAAGGGAATGTCATACACGATGTGTTCAGGAAAGTTTTCAATTGACAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTGAAAGTCAAGTATGAAGGTGCTGGAGCTCCGTGTAAAGTCCCCATAGAGATAAGAGATGTAAACAAGGAAAAAGTGGTTGGGCGTATCATCTCATCCACCCCTTTGGCTGAGAATACCAACAGTGTAACCAACATAGAATTAGAACCCCCCTTTGGGGACAGCTACATAGTGATAGGTGTTGGAAACAGCGCATTAACACTCCATTGGTTCAGGAAAGGGAGTTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGTGCAAAACGAATGGCCATTCTAGGTGAAACAGCTTGGGATTTTGGTTCCGTTGGTGGACTGTTCACATCATTGGGAAAGGCTGTGCACCAGGTTTTTGGAAGTGTGTATACAACCATGTTTGGAGGAGTCTCATGGATGATTAGAATCCTAATTGGGTTCTTAGTGTTGTGGATTGGCACGAACTCAAGGAACACTTCAATGGCTATGACGTGCATAGCTGTTGGAGGAATCACTCTGTTTCTGGGCTTCACAGTTCAAGCA",
           "ATGCGATGCGTGGGAGTGGGGAACAGAGACTTTGTGGAAGGAGTCTCAGGTGGAGCATGGGTCGATTTGGTGCTAGAACATGGAGGATGTGTCACAACCATGGCCCAGGGAAAACCAACCTTGGATTTTGAACTGATCAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCGATATCAAACATAACCACGGCAACAAGATGTCCAACGCAAGGAGAACCTTATCTCAAAGAGGAACAAGATCAACAGTACATCTGCCGGAGAGATGTGGTAGACAGAGGGTGGGGCAATGGCTGTGGCTTGTTTGGGAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGCTCGGGGAAGATAACAGGCAATTTGGTCCAAATTGAGAACCTTGAATACACAGTGGTTGTAACAGTCCACAATGGAGACACCCATGCAGTAGGAAATGATACATCCAACCATGGAGTGACAGCCACGATAACCCCCAGGTCACCATCGGTAGAAGTTAAATTACCGGATTATGGAGAATTAACACTCGATTGTGAACCCAGGTCCGGAATTGATTTTAATGAGATGATTCTGATGAAAATGAAAAAGAAAACGTGGCTTGTGCACAAGCAATGGTTTTTGGATCTACCTCTACCATGGGCAGCAGGAGCAGATACATCAGAAGTTCATTGGAATTACAAAGAGAGAATGGTGACATTCAAGGTTCCTCATGCCAAGAGACAGGATGTGACAGTGCTAGGATCTCAGGAAGGAGCCATGCATTCTGCCCTCACCGGAGCTACAGAAGTGGATTCCGGTGATGGAAACCACATGTTTGCAGGACATCTGAAATGCAAAGTTCGCATGGAGAAATTGAGAATTAAGGGAATGTCATACACGATGTGCTCAGGAAAGTTCTCAATTGACAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTAAAAGTCAAATATGAGGGTGCTGGAGCTCCATGTAAAGTTCCCATAGAGATAAGAGATGTGAACAAGGAAAAAGTGGTAGGGCGTATCATCTCATCTACCCCTTTGGCTGAGAACACCAACAGTGTAACCAACATAGAATTAGAACCCCCTTTTGGGGACAGCTACATAGTAATAGGTGTTGGAGACAGTGCATTAACACTCCATTGGTTCAGGAAAGGGAGTTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGTGCAAAGCGAATGGCCATTCTAGGTGAAACAGCCTGGGATTTTGGTTCGGTTGGTGGACTGCTCACATCATTGGGAAAGGCTGTACACCAGGTTTTTGGTAGTGTGTATACAACTATGTTTGGAGGAGTCTCATGGATGGTTAGAATCCTAATTGGGTTCTTAGTGTTGTGGATTGGCACGAATTCGAGAAACACCTCAATGGCAATGACGTGCATAGCTGTTGGAGGAATCACTCTGTTTCTGGGTTTCACAGTTCACGCA",
           "ATGCGATGCGTGGGAGTGGGGAACAGAGACTTTGTGGAAGGAGTCTCAGGTGGAGCATGGGTCGATTTGGTGCTAGAACATGGAGGATGTGTCACAACCATGGCCCAGGGAAAACCAACCTTGGATTTTGAACTGATCAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCGATATCAAACATAACCACGGCAACAAGATGTCCAACGCAAGGAGAACCTTATCTCAAAGAGGAACAAGATCAACAGTACATCTGCCGGAGAGATGTGGTAGACAGAGGGTGGGGCAATGGCTGTGGCTTGCTTGGGAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGCTCGGGGAAGATAACAGGCAATTTGGTCCAAATTGAGAACCTTGAATACACAGTAGTTGTAACAGTCCACAATGGAGACACCCATGCAGTAGGAAATGACATATCCAACCATGGAGTGACAGCCACGATAACCCCCAGGTCACCATCGGTAGAAGTTAAATTACCGGATTATGGAGAATTAACACTCGATTGTGAACCCAGGTCCGGAATTGATTTTAATGAGATGATTCTGATGAAAATGAAAAAGAAAACGTGGCTTGTGCACAAGCAATGGTTTTTGGATCTACCTCTACCATGGGCAGCAGGAGCAGACACATCAGAAGTTCATTGGAATTACAAAGAGAGAATGGTGACATTCAAGGTTCCTCATGCCAAGAGACAGGATGTGACAGTGCTAGGATCTCAGGAAGGAGCCATGCATTCTGCCCTCACCGGAGCTACAGAAGTGGATTCCGGTGATGGAAACCACATGTATGCAGGACATCTGAAATGCAAAGTTCGCATGGAGAAATTGAGAATTAAGGGAATGTCATACACGATGTGCTCAGGAAAGTTCTCAATTGACAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTAAAAGTCAAGTATGAGGGTGCTGGAGCTCCATGTAAAGTTCCCATAGAGATAAGAGATGTGAACAAGGAAAAAGTGGTAGGGCGCATCATCTCATCTACCCCTTTGGCTGAGTATACCAACAGTGTAACCAACATAGAATTAGAACCCCCCTTTGGGGACAGCTACATAGTAATAGGTGTTGGAGACAGTGCATTAACACTCCATTGGTTCAGGAAAGGGAGTTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGCGCAAAGCGAATGGCCATTCTAGGTGAAACAGCCTGGGATTTTGGTTCTGTTGGTGGACTGCTCACATCATTGGGAAAGGCTGTACACCAGGTTTTTGGTAGTGTGTATACAACTATGTTTGGAGGAGTCTCATGGATGGTTAGAATCCTAATTGGGTTCTTAGTGTTGTGGATTGGCACGAATTCGAGAAACACCTCAATGGCAATGACGTGCATAGCTGTTGGAGGAATCACTCTGTTTCTGGGTTTCACAGTTCACGCA",
           "ATGCGATGCGTAGGAGTGGGGAACAGAGACTTTGTGGAAGGAGTCTCAGGTGGAGCATGGGTCGACTTAGTGCTAGAACATGGAGGATGTGTCACAACCATGGCCCAAGGAAAACCAACCTTGGATTTTGAACTGATCAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCGATATCAAACATAACCACGGCAACAAGATGCCCAACGCAAGGAGAACCTTATCTCAAAGAGGAACAAGATCAACAGTACATTTGCCGGAGAGATGTGGTAGACAGAGGGTGGGGCAATGGCTGTGGCTTGTTTGGGAAAGGAGGAGTTGTGACATGTGCGAAGTTCTCATGCTCGGGAAAGATAACAGGCAATTTGGTCCAAATTGAGAACCTTGAATATACAGTGGTTGTAACAGTCCACAATGGAGACACCCATGCAGTAGGAAATGACACATCCAACCATGGAGTGACAGCCACGATAACCCCTAGGTCACCATCGGTAGAAGTTAAATTACCGGATTATGGAGAATTAACACTTGATTGTGAACCCAGGTCCGGAATTGACTTTAATGAGATGATTCTGATGAAAATGAAAAAGAAAACGTGGCTTGTGCACAAGCAATGGTTTCTGGATCTGCCTCTACCATGGGCAGCAGGAGCAGATACATCAGAAGTTCATTGGAATTACAAAGAGAGAATGGTGACATTCAAGGTTCCTCATGCCAAGAGACAGGATGTGACAGTGCTAGGATCTCAGGAAGGAGCCATGCATTCTGCCCTCACCGGAGCTACAGAAGTGGATTCCGGTGATGGAAACCACATGTTTGCAGGACATCTGAAATGCAAAGTTCGCATGGAGAAATTGAGAATCAAGGGAATGTCATACACGATGTGCTCAGGGAAGTTCTCAATTGACAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTTAAAGTCAAATATGAAGGTGCTGGAGCTCCGTGCAAAGTTCCCATAGAGATAAGAGATGTGAACAGGGAAAAAGTGGTAGGGCGTGTCATCTCATCTACCCCTTTGGCCGAGAATACCAACAGTGTAACCAACATAGAATTAGAACCCCCTTTTGGGGACAGCTACATAGTAATAGGTGTTGGAGACAGTGCATTAACACTCCATTGGTTCAGGAAAGGGAGTTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGTGCAAAGCGAATGGCCATTCTAGGTGAAACAGCCTGGGATTTTGGTTCTGTTGGTGGACTGCTCACATCATTGGGAAAGGCTGTACACCAGGTTTTTGGTAGTGTGTATACAACTATGTTTGGAGGAGTCTCATGGATGGTTAGAATCCTAATTGGGTTCCTAGTGTTGTGGATTGGCACGAATTCGAGAAACACCTCAATGGCAATGACGTGCATAGCTGTTGGAGGAATCACTCTGTTCATGGGTTTCACAGTTCACGCA",
           "ATGCGATGCGTAGGAGTAGGAAACAGAGACTTTGTGGAAGGAGTCTCAGGTGGAGCATGGGTCGACCTGGTGCTAGAACATGGAGGATGCGTCACAACCATGGCCCAGGGAAAACCAACCTTGGATTTTGAACTGACTAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCAATATCAAACATAACTACGGCAACAAGATGTCCAACGCAAGGGGAGCCCTATCTGAAAGAGGAACAGGACCAACAGTACATTTGCCGGAGAGATGTGGTAGACAGAGGGTGGGGCAATGGCTGTGGCTTGTTTGGAAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGTTCGGGGAAGATAACAGGCAATCTGGTCCAAATTGAGAACCTTGAGTACACAGTGGTTGTAACAGTCCACAATGGAGACACCCATGCAGTAGGAAATGACACATCCAATCATGGAGTTACAGCCACGATAACTCCCAGGTCACCATCGGTGGAAGTCAAATTGCCGGACTATGGAGAACTAACACTCGATTGTGAACCCAGGTCTGGAATTGACTTTAATGAGATGATTCTAATGAAAATGAAAAAGAAAACATGGCTCGTGCATAAGCAATGGTTTTTGGATCTGCCTCTTCCATGGACAGCAGGAGCAGACACATCAGAGGTTCACTGGAATTACAAAGAGAGAATGGTGACATTTAAGGTTCCTCATGCCAAGAGACAGGATGTGACAGTGCTGGGATCTCAGGAAGGAGCCATGCATTCTGCCCTCGCTGGAGCCACAGAAGTGGACTCCGGTGATGGAAATCACATGTTTGCAGGACATCTCAAGTGCAAAGTCCGTATGGAGAAATTGAGAATCAAGGGAATGTCATACACGATGTGTTCAGGAAAGTTTTCAATCGACAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTGAAAGTCAAGTATGAAGGTGCTGGAGCTCCGTGTAAAGTCCCCATAGAGATAAGAGATGTAAACAAGGAAAAAGTGGTTGGGCGCGTCATCTCATCCACCCCTTTGGCTGAGAATACCAACAGTGTAACCAACATAGAATTAGAACCCCCCTTTGGGGACAGCTACATAGTGATAGGTGTTGGAAACAGCGCATTAACACTCCATTGGTTCAGGAAAGGGAGTTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGTGCAAAACGAATGGCCATTCTAGGTGAAACAGCTTGGGATTTTGGTTCCGTTGGTGGACTGTTCACATCATTGGGAAAGGCTGTGCACCAGGTTTTTGGAAGCGTGTATACAACCATGTTTGGAGGAGTCTCATGGATGATTAGAATCCTAATTGGGTTCTTAGTGTTGTGGATTGGCACGAACTCAAGGAACACTTCAATGGCTATGACGTGCATAGCTGTTGGAGGTATCACTCTGTTTCTGGGCTTCACAGTTCAAGCG",
           "ATGCGATGCGTGGGAGTGGGGAACAGAGACTTTGTGGAAGGAGTCTCAGGTGGAGCATGGGTCGATCTGGTGCTAGAACATGGAGGATGTGTCACAACCATGGCCCAAGGAAAACCAACCTTGGATTTTGAACTGATCAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCCATATCAAACATAACCACGGCAACAAGATGTCCAACGCAAGGAGAGCCCTGTCTCAAAGAGGAACAGGATCAACAGTACATCTGCCGGAGAGACGTGGTAGACAGAGGGTGGGGCAATGGCTGTGGCCTGCTTGGGAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGCTCGGGGAAGATAACAGGCAACTTAGTCCGAATTGAGAACCTTGAATACACAGTGGTTGTAACAGTCCACAATGGAGACACCCATGCAGTAGGAAATGACACATCCAACCACGGAGTGACAGCCACGATAACCCCCAGGTCACCATCGGTAGAAGTTAAATTACCGGACTATGGAGAATTGACACTCGATTGTGAACCCAGGTCCGGAATTGATTTTAATGAGATGATTCTGATGGAAATGAGAAAGAAAACGTGGCTTGTGCACAAGCAATGGTTTTTGGATCTACCTCTACCATGGACAGCAGGAGCAGACACGTCAGAAGTTCATTGGAATCACAAAGAGAGAATGGTGACGTTCAAGGTCCCTCATGCCAAGAGACAGGATGTGACAGTGCTAGGATCTCAGGAAGGAGCCATGCATTCAGCCCTCACCGGAGCCACAGAAGTGGATTCCGGTGATGGAAACCACATGTTTGCAGGACATTTGAAGTGCAAAGTTCGCATGGAGAAATTGAGAATCAAGGGAATGTCATACACGATGTGCTCAGGAAAGTTCTCAATTGATAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTAAAAGTCAAGTATGAGGGTGCCGGAGCTCCATGTAAAGTTCCCATAGAGATAAGAGATGTGAACAAGGAAAAAGTGGTTGGGCGCATCATCTCATCTACCCCTTTTGCTGAGAATACCAACAGTGTGACCAACATAGAATTGGAACCCCCCTTTGGGGATAGCTACATAGTAATAGGTGTAGGAAACAGTGCATTAACACTCCATTGGTTTAGGAAAGGGAGTTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGCGCAAAGCGCATGGCCATTCTAGGTGAAACAGCTTGGGATTTTGGTTCTGTTGGTGGACTGCTCACATCATTGGGAAAGGCTGTACACCAGGTTTTTGGTAGTGTGTATACAACTATGTTTGGAGGAGTCTCATGGATGGTTAGAATCCTAATCGGGTTCTTAGTATTGTGGATTGGCACGAATTCAAGAAACACTTCAATGGCAATGTCGTGCATAGCTGTTGGAGGAATTACTTTGTTTCTGGGTTTCACAGTTCATGCA",
           "ATGCGATGCGTAGGAGTAGGAAACAGAGACTTTGTGGAAGGAGTCTCAGGTGGAGCATGGGTCGATCTGGTGCTAGAACATGGAGGATGCGTCACAACCATGGCCCAGGGAAAACCAACCTTGGATTTTGAACTGACTAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCAATATCAAACATAACTACGGCAACAAGATGTCCAACGCAAGGAGAGCCTTATCTGAAAGAGGAACAGGACCAACAGTACATTTGCCGGAGAGATGTGGTAGACAGAGGGTGGGGCAATGGCTGTGGCTTGTTTGGAAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGTTCGGGGAAGATAACAGGCAATTTGGTCCAAATTGAGAACCTTGAATACACAGTGGTTGTAACAGTCCACAATGGAGACACCCATGCAGTAGGAAATGACACATCCAATCATGGAGTTACAGCCACGATAACTCCCAGGTCACCATCGGTGGAAGTCAAATTGCCGGACTATGGAGAACTAACACTCGATTGTGAACCCAGGTCTGGAATTGACTTTAATGAGATGATCCTGATGAAAATGAGAAAGAAGACATGGCTCGTGCATAAGCAATGGTTTTTGGATCTGCCTCTTCCATGGACAGCAGGAGCAGACACATCAGAGGTTCACTGGAATTACAAAGAGAGAATGGTGACATTTAAGGTTCCTCATGCCAAGAGACAGGATGTGACAGTGCTGGGATCTCAGGAAGGAGCCATGCATTCTGCCCTCGCTGGAGCCACAGAAGTGGACTCCGGTGATGGAAATCACATGTTTGCAGGACATCTCAAGTGCAAAGTCCGTATGGAGAAATTGAGAATCAAGGGAATGTCATACACGATGTGTTCAGGAAAGTTTTCAATTGACAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTGAAAGTCAAGTATGAAGGTGCTGGAGCTCCGTGTAAAGTCCCCATAGAGATAAGAGATGTAAACAAGGAAAAAGTGGTTGGGCGTATCATCTCATCCACCCCTTTGGCTGAGAATACCAACAGTGTAACCAACATAGAATTAGAACCCCCCTTTGGGGACAGCTACATAGTGATAGGTGTTGGAAACAGCGCATTAACACTCCATTGGTTCAGGAAAGGGAGTTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGTGCAAAACGAATGGCCATTCTAGGTGAAACAGCTTGGGATTTTGGTTCCGTTGGTGGACTGTTCACATCATTGGGAAAGGCTGTGCACCAGGTTTTTGGAAGTGTGTATACAACCATGTTTGGAGGAGTCTCATGGATGATTAGAATCCTAATTGGGTTCTTAGTGTTGTGGATTGGCACGAACTCAAGGAACACTTCAATGGCTATGACGTGCATAGCTGTTGGAGGAATCACTCTGTTTCTGGGCTTCACAGTTCAAGCA",
           "ATGCGATGCGTAGGAGTAGGAAACAGAGACTTTGTGGAAGGAGTTTCAGGTGGAGCATGGGTCGATTTGGTGCTAGAACATGGAGGATGCGTCACAACCATGGCCCAGGGAAAACCAACCTTGGATTTTGAACTGACTAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCAATATCAAACATAACTACGGCAACAAGATGTCCAACGCAAGGAGAGCCTTATCTGAAAGAGGAACAGGACCAACAGTACATTTGCCGGAGAGATGTGGTAGACAGAGGGTGGGGCAATGGCTGTGGCTTGTTTGGAAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGTTCGGGGAAGATAACAGGCAATCTGGTCCAAATTGAGAACCTTGAATACACAGTGGTCATAACAGTCCACAATGGAGACACCCATGCAGTAGGAAATGACACATCCAATCATGGAGTTACAGCCACGATAACTCCCAGGTCACCATCGGTGGAAGTCAAATTGCCGGACTATGGAGAACTAACACTCGATTGTGAACCCAGGTCTGGAATTGACTTTAATGAGATGATTCTGATGAAAATGAAAAAGAAAACATGGCTCGTGCATAAGCAATGGTTTTTGGATCTGCCTCTTCCATGGACAGCAGGAGCAGACACAACAGAGGTTCACTGGAATTACAAAGAGAGAATGGTGACATTTAAGGTTCCTCATGCCAAGAGACAGGATGTGACAGTACTGGGATCTCAGGAAGGAGCCATGCATTCTGCCCTAGCTGGAGCTACAGAAGTGGACTCCGGTGATGGGAATCACATGTTTGCAGGACATCTCAAGTGCAAAATCCGTATGGAGAAATTGAGAATCAAGGGAATGTCATACACGATGTGTTCAGGAAAGTTTTCAATTGACAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTGAAAGTCAAGTATGAAGGTGCTGGAGCTCCGTGTAAAGTCCCCATAGAGATAAGAGATGTAAACAAGGAAAAAGTGGTTGGGCGTGTCATCTCATCCACCCCTTTGGCTGAGAATACCAACAGTGTAACCAACATAGAATTAGAACCCCCCTTTGGGGACAGCTACATAGTGATAGGTGTTGGAAACAGCGCATTAACACTCCATTGGTTCAGGAAAGGGAGTTCCATCGGCAAGATGTTTGAGTCCACATACAGAGGTGCAAAACGAATGGCCATTCTAGGTGAAACAGCTTGGGATTTTGGTTCCGTTGGTGGACTGTTCACATCATTGGGAAAGGCTGTGCACCAGGTTTTTGGAAGTGTGTATACAACCATGTTTGGAGGAGTCTCATGGATGATTAGAATCCTAATTGGGTTCTTAGTGTTGTGGATTGGCACGAACTCAAGGAACACTTCAATGGCTATGACGTGCATAGCTGTTGGAGGAATCACTCTGTTTCTGGGCTTCACAGTTCAAGCA",
           "ATGCGATGCGTAGGAGTGGGGAACAGGGACTTTGTGGAAGGAGTCTCAGGTGGAGCATGGGTCGATCTGGTGCTAGAACATGGAGGATGTGTCACAACCATGGCTCAAGGAAAACCAACCTTGGATTTTGAACTGATCAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCGATATCAAACATAACCACGGCGACAAGATGTCCAACGCAAGGAGAGCCTTATCTCAAAGAGGAACAAGATCAACAGTACATCTGCCGGAGAGATGTGGTAGACAGAGGGTGGGGCAATGGCTGTGGCTTACTTGGAAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGCTCGGGGAAGATAACAGGCAACTTGGTTCGAATTGAGAACCTTGAATACACAGTGGTTGTGACAGTCCACAACGGAGACACCCATGCAGTAGGAAATGACATATCCAACCATGGAGTGACAGCCACGATAACTCCCAGGTCACCATCGGTAGAAGTCAAATTACCGGATTATGGAGAATTAACGCTCGATTGTGAACCCAGGTCCGGAATTGATTTTAATGAGATGATTCTGATGGAAATGAGAAAGAAGACGTGGCTTGTGCACAAGCAATGGTTTTTGGATCTACCTCTACCATGGACAGCAGGAGCAGACACAGCAGAAGTTCATTGGAATTACAAAGAGAGAATGGTGACATTCAAGGTCCCTCATGCTAAGAGACAAGATGTGACAGTGCTAGGATCTCAGGAAGGAGCCATGCATTCTGCCCTCACCGGAGCTACAGAAGTGGATTCCGGTGATGGAAACCACATGTTTGCAGGACATCTTAAGTGCAAGGTTCGTATGGAGAAATTGAGAATCAAGGGAATGTCATACACGATGTGCTCAGGAAAGTTCTCAATTGACAAAGAGATGGCAGAAACACAGCATGGGACAACAGTAGTGAAAGTTAAGTATGAAGGCGCTGAAGCTCCATGTAAAATCCCCATAGAGATAAGAGATGTGAACAAGGAAAAAGTTGTTGGGCGCATCATCTCATCCACTCCTTTGGCTGAAAACACCAACAGCGTGACCAATATAGAATTAGAACCTCCCTTTGGGGACAGCTACATAGTAATAGGTGTTGGAGACAGTGCATTAACACTCCATTGGTTCAGGAAGGGGAGTTCCATTGGCAAGATGTTTGAGTCTACATATAGAGGAGCAAAGCGAATGGCCATTCTAGGTGAAACAGCTTGGGACTTTGGCTCTGTTGGCGGACTGTCTACATCATTGGGAAAGGCTGTACACCAGGTTTTTGGTAGTGTGTACACAACCATGTTTGGAGGAGTCTCATGGATGGTTAGAATCCTAATCGGGCTCTTGGTGTTGTGGATTGGCACAAATTCAAGAAACACCTCAATGGCAATGACGTGCATAGCTGTTGGAGGAATCACTCTATTTCTGGGTTTCACAGCTCACGCA",
           "ATGCGATGCGTAGGAGTGGGGAACAGAGACTTTGTAGAAGGAGTCTCAGGTGGAGCATGGGTCGATCTGGTGCTAGAACATGGAGGATGTGTCACAACCATGGCCCAGGGAAAACCAACCTTGGATTTTGAACTGATCAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCCATATCAAACATAACCACGGCAACAAGATGTCCAACGCAAGGAGAGCCTTACCTCAAAGAGGAACAAGATCAACAGTACATCTGCCGGAGAGACGTGGTAGATAGAGGGTGGGGCAACGGCTGTGGCTTGCTTGGGAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGCTCGGGGAAGATAACAGGCAACTTAGTCCAAATTGAGAACCTTGAATACACAGTGGTTGTAACAGTCCACAATGGAGACACCCATGCTGTAGGAAATGATACATCCAACCACGGAGTGACAGCCACGATAACCCCCAGGTCACCATCGGTAGAAGTTAAATTACCGGACTATGGAGAATTAACACTTGATTGTGAACCTAGGTCCGGAATTGACTTTAATGAGATGATTCTGATGAAAATGAAAAAGAAAACGTGGCTCGTGCACAAGCAATGGTTTTTGGATCTACCTCTACCATGGACAGCAGGAGCAGACACGTCAGAAGTTCACTGGAATCACAAAGAGAGAATGGTGACATTCAAGGTTCCTCATGCCAAGAGACAGGATGTGACAGTGCTAGGATCTCAGGAAGGAGCTATGCATTCAGCCCTCACCGGAGCCACAGAAGTGGATTCCGGTGATGGAAACCATATGTTTGCAGGACATCTTAAGTGTAAAGTTCGTATGGAGAAATTGAGGATCAAGGGAATGTCATACACGATGTGCTCAGGAAAGTTCTCAATTGATAAAGAGATGGCAGAAACACAGCATGGGACAACAGTGGTAAAAGTCAAGTATGAAGGTGCTGGAGCTCCATGTAAAGTCCCCATAGAGATAAGAGATGTGAACAAGGAAAAAGTGGTTGGGCGTATCATCTCATCTACCCCTTTTGCTGAGAATACCAACAGTGTGACCAATATAGAATTGGAACCCCCTTTTGGGGATAGCTACATAGTAATAGGTGTAGGAGACAGTGCATTAACACTCCATTGGTTCAGGAAAGGGAGCTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGCGCAAAGCGCATGGCCATTCTAGGTGAAACAGCTTGGGATTTTGGTTCTGTCGGTGGACTGCTCACATCATTGGGAAAGGCTGTACACCAGGTTTTTGGCAGTGTGTATACAACTATGTTTGGAGGAGTCTCATGGATGGTTAGAATCCTAATCGGGTTCTTAGTGTTGTGGATTGGCACGAATTCAAGAAACACTTCAATGGCAATGTCGTGCATAGCTGTTGGAGGAATCACTCTGTTTCTGGGCTTCACAGTTCATGCA",
           "ATGCGATGCGTAGGAGTAGGGAACAGAGACTTTGTAGAAGGAGTCTCAGGTGGAGCATGGGTCGATCTGGTGCTAGAACATGGAGGATGTGTCACAACCATGGCCCAGGGAAAACCAACCCTGGATTTTGAACTGATCAAGACAACAGCCAAGGAAGTGGCTCTGTTAAGAACCTATTGCATTGAAGCCTCCATATCAAACATAACCACGGCAACAAGATGTCCAACGCAAGGAGAGCCTTACCTCAAAGAGGAACAAGATCAACAGTACATCTGCCGGAGAGACGTGGTAGACAGAGGGTGGGGCAACGGCTGTGGCTTGTTTGGGAAAGGAGGAGTTGTGACATGTGCGAAGTTTTCATGCTCGGGGAAGATAACGGGCAACTTAGTCCAAATTGAGAACCTTGAATACACAGTGGTTGTAACAGTCCACAATGGAGACACCCATGCTGTAGGAAATGATACATCCAACCACGGAGTGACAGCCACGATAACCCCCAGGTCACCATCGGTAGAAGTTAAATTACCGGACTATGGAGAATTAACACTTGATTGTGAACCTAGGTCCGGAATTGATTTTAATGAGATGATTCTGATGAAAATGAAAAAGAAAACGTGGCTCGTGCACAAGCAATGGTTTTTGGATCTACCTCTACCATGGACAGCAGGAGCAGACACGTCAGAAGTTCACTGGAATCACAAAGAGAGAATGGTGACATTCAAGGTTCCTCATGCCAAGAGACAGGATGTGACAGTGCTAGGATCTCAGGAAGGAGCCATGCATTCAGCCCTCGCCGGAGCCACAGAAGTGGATTCCGGTGATGGAAACCATATGTTTGCAGGACACTTAAAGTGTAAAGTTCGTATGGAGAAATTGAGAATCAAGGGAATGTCATACACGATGTGCTCAGGAAAGTTCTCAATTGATAAAGAGATGGCAGAAACACAGCATGGAACAACAGTGGTAAAGGTCAAGTATGAAGGCACTGGAGCTCCATGTAAAGTCCCCATAGAGATAAGAGATGTGAACAAGGAAAAAGTGGTTGGGCGTATCATCTCATCTACCCCTTTGGCTGAGAATACCAACAGTGTGACCAATATAGAATTGGAACCTCCTTTTGGGGATAGCTACATAGTGATAGGTGTGGGAGACAGTGCATTAACACTCCATTGGTTCAGGAAAGGGAGCTCCATTGGCAAGATGTTTGAGTCCACATACAGAGGTGCAAAGCGCATGGCCATTCTAGGTGAAACAGCTTGGGATTTTGGTTCTGTTGGCGGACTGCTCACATCATTGGGAAAGGCTGTACACCAGGTTTTTGGCAGTGTGTATACAACTATGTTTGGAGGAGTCTCATGGATGGTTAGAATCCTAATCGGGTTCTTAGTGTTGTGGATTGGCACGAATTCAAGAAACACTTCAATGGCTATGTCGTGCATAGCTGTTGGAGGAATTACTCTGTTTCTGGGCTTCACAGTTCATGCA"},
          {1982.0, 1983.0, 1994.0, 1976.0, 1977.0, 1984.0, 1981.0, 1964.0, 1956.0, 1984.0, 1986.0, 1978.0, 1979.0, 1985.0, 1963.0, 1978.0, 1984.0},
          {false}//forwards
          };


    protected ConstantPopulationModel createRandomInitialTree(Parameter popSize) {        
        ConstantPopulationModel startingTree = new ConstantPopulationModel(popSize, Units.Type.YEARS);
        ConstantPopulation constant = (ConstantPopulation) startingTree.getDemographicFunction();

        createTreeModel(constant);

        return startingTree;
    }

    private void createTreeModel (ConstantPopulation constant) {
        CoalescentSimulator simulator = new CoalescentSimulator();
        Tree tree = simulator.simulateTree(alignment, constant);
        treeModel = new TreeModel(tree);//treeModel
    }

}

