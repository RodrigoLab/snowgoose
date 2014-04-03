package test.srp.spectrum.treelikelihood;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.apache.commons.math3.stat.StatUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.MCMCSetupHelper;
import srp.core.MCMCSetupHelperSpectrum;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.HaplotypeModelUtils;
import srp.operator.spectrum.DeltaExchangeSingleSpectrumOperator;
import srp.shortreads.AlignmentMapping;
import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumAlignmentUtils;
import srp.spectrum.treelikelihood.SpectrumTreeLikelihood;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.ext.TreeLikelihoodExt;
import dr.inference.markovchain.MarkovChain;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.math.MathUtils;

public class SpectrumTreeLikelihoodTest {

	private static final double THRESHOLD = MarkovChain.EVALUATION_TEST_THRESHOLD;
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
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
	public void testSetting() throws Exception {
		// Sub model
		Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25,
				0.25, 0.25 });
		Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100);

		FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
		HKY hky = new HKY(kappa, f);

		// siteModel
		GammaSiteModel siteModel = new GammaSiteModel(hky);
		Parameter mu = new Parameter.Default(
				GammaSiteModelParser.MUTATION_RATE, 1.0, 0,
				Double.POSITIVE_INFINITY);
		siteModel.setMutationRateParameter(mu);

		// treeLikelihood
		Taxon[] taxa = new Taxon[4];
		for (int i = 0; i < taxa.length; i++) {
			taxa[i] = new Taxon("taxa_"+i);
		}
		SimpleAlignment alignment = new SimpleAlignment();
		alignment.addSequence(new Sequence(taxa[0], "A"));
		alignment.addSequence(new Sequence(taxa[1], "A"));
		alignment.addSequence(new Sequence(taxa[2], "C"));
		alignment.addSequence(new Sequence(taxa[3], "C"));
		SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
		TreeModel treeModel = createTreeModel(taxa, 0.01);
		TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel,
				siteModel, null, null, false, false, true, false, false);
		double likelihood = treeLikelihood.getLogLikelihood();

		
		String[] seqs = new String[] {
				"A",
				"A",
				"C",
				"C"};
			
		alignment = AlignmentUtils.createAlignment(seqs);
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(alignment);
		SpectrumTreeLikelihood spectrumTreeLikelihood = new SpectrumTreeLikelihood(spectrumModel, treeModel,
				siteModel, null, false, false, true, false, false);
		double srpTreeLikelihood = spectrumTreeLikelihood.getLogLikelihood();
		assertEquals(likelihood, srpTreeLikelihood, 0);
	}

	@Test
	public void testSpectrumModel() throws Exception {
		// Sub model
		Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25,
				0.25, 0.25 });
		Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100);

		FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
		HKY hky = new HKY(kappa, f);

		// siteModel
		GammaSiteModel siteModel = new GammaSiteModel(hky);
		Parameter mu = new Parameter.Default(
				GammaSiteModelParser.MUTATION_RATE, 1.0, 0,
				Double.POSITIVE_INFINITY);
		siteModel.setMutationRateParameter(mu);

		// treeLikelihood
		int taxaCount = 2;
		Taxon[] taxa = new Taxon[taxaCount];
		for (int i = 0; i < taxa.length; i++) {
			taxa[i] = new Taxon("taxa_"+i);
		}
		TreeModel treeModel = createTreeModel2Taxa(taxa, 0.01);

		String[] seqs = new String[] {
				"A",
				"A",
				"C",
				"C"};
			
		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, 1);
		spectrumModel.removeSpectrum(0);
					
//		4	[0.9802314196914651, 1.096410868086927E-5, 1.096410868086935E-5, 1.096410868086927E-5]
//		[1.0, 0.0, 0.0, 0.0]	[1.0, 0.0, 0.0, 0.0]
//		5	[1.0964108663326694E-5, 0.9802314196862197, 1.0964108663326694E-5, 1.0964108663326615E-5]
//		[0.0, 1.0, 0.0, 0.0]	[0.0, 1.0, 0.0, 0.0]
//		6	[1.2757307680901868E-5, 1.2757307699647983E-5, 1.6938161666098327E-10, 1.6938161666108597E-10]
//		[0.9802314196914651, 1.096410868086927E-5, 1.096410868086935E-5, 1.096410868086927E-5]	[1.0964108663326694E-5, 0.9802314196862197, 1.0964108663326694E-5, 1.0964108663326615E-5]

//		[0.0031605887470258776, 0.0031605887495542917, 1.0605960978378439E-5, 1.0605960978378359E-5]
//				[0.9802314196914651, 1.096410868086927E-5, 1.096410868086935E-5, 1.096410868086927E-5]	
//						[1.0964108663326694E-5, 0.9802314196862197, 1.0964108663326694E-5, 1.0964108663326615E-5]

		int stateCount = 4;
		
		double[] node1 = new double[]{0.9802314196914651, 1.096410868086927E-5, 1.096410868086935E-5, 1.096410868086927E-5};
		Spectrum s = new Spectrum(spectrumLength);
		for (int i = 0; i < s.getLength(); i++) {
			SpectraParameter spectra = s.getSpectra(i); 
			for (int j = 0; j < stateCount; j++) {
				spectra.setFrequency(j, node1[j]);
			}
		}
		s.setTaxon(taxa[0]);
		spectrumModel.addSpectrum(s);
		
		s = new Spectrum(spectrumLength);
		double[] node2 = new double[]{1.0964108663326694E-5, 0.9802314196862197, 1.0964108663326694E-5, 1.0964108663326615E-5};
		for (int i = 0; i < s.getLength(); i++) {
			SpectraParameter spectra = s.getSpectra(i); 
			for (int j = 0; j < stateCount; j++) {
				spectra.setFrequency(j, node2[j]);
			}
		}
		s.setTaxon(taxa[1]);
		spectrumModel.addSpectrum(s);
		
		
		SpectrumTreeLikelihood treeLikelihood = new SpectrumTreeLikelihood(spectrumModel, treeModel,
				siteModel, null, false, false, true, false, false);

		System.out.println(treeLikelihood.getLogLikelihood());
		
		double[] expecteds = new double[]{0.0031605887470258776, 0.0031605887495542917, 1.0605960978378439E-5, 1.0605960978378359E-5};
		assertArrayEquals(expecteds, treeLikelihood.getRootPartials(), 0);
		
		double[] rootPartial = treeLikelihood.getRootPartials();
		double expected = -6.446794062758868;
		double logLikelihood = Math.log(StatUtils.sum(rootPartial)*0.25);
		assertEquals(expected, logLikelihood, 0) ;
//		System.out.println(Arrays.toString(rootPartial));
//		System.out.println(StatUtils.sum(rootPartial)*0.25 +"\t"+ Math.log(StatUtils.sum(rootPartial)*0.25));
//		System.out.println(Arrays.toString(treeLikelihood.getPatternLogLikelihoods()));
//		LikelihoodCore likelihoodCore = treeLikelihood.getLikelihoodCore();
//		likelihoodCore.calculatePartials(nodeIndex1, nodeIndex2, nodeIndex3)
	}
	
	@Test
	public void testUpdateSpectrum() {

		SiteModel siteModel = MCMCSetupHelper.setupSiteModel();
		int taxaCount = 8;
		Taxon[] taxa = new Taxon[taxaCount];
		for (int i = 0; i < taxa.length; i++) {
			taxa[i] = new Taxon("taxa_"+i);
		}

		TreeModel treeModel = createTreeModel(taxa, 0.01);
		
		String[] seqs = new String[] {"ACGTTTA"};
			
		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, taxaCount);
		SpectrumTreeLikelihood spectrumTreeLikelihood = new SpectrumTreeLikelihood(spectrumModel, treeModel,
				siteModel, null, false, false, true, false, false);

		MCMCOperator op = new DeltaExchangeSingleSpectrumOperator(spectrumModel, 0.1, null);
		double likelihood = spectrumTreeLikelihood.getLogLikelihood();
		double expected = 0;
		for (int i = 0; i < 1e4; i++) {
			try {
				spectrumTreeLikelihood.storeModelState();
				
				op.operate();
				SpectrumAlignmentModel newSpectrumModel = SpectrumAlignmentModel.duplicateSpectrumAlignmentModel(spectrumModel);
				SpectrumTreeLikelihood newSpectrumTreeLikelihood = new SpectrumTreeLikelihood(newSpectrumModel, treeModel,
						siteModel, null, false, false, true, false, true);
				expected = newSpectrumTreeLikelihood.getLogLikelihood();
				likelihood = spectrumTreeLikelihood.getLogLikelihood();
				assertEquals(expected, likelihood, 1e-10);
				double rand = MathUtils.nextDouble();

				if(rand>0.5){
					op.accept(0);
					spectrumTreeLikelihood.acceptModelState();
				}
				else{
					op.reject();
					spectrumTreeLikelihood.restoreModelState();
				}

			} catch (OperatorFailedException e){
				op.reject();
			} catch (AssertionError e){
				System.out.println(i);
				System.out.println(likelihood);
				System.out.println(expected);
				SpectrumTreeLikelihood newSpectrumTreeLikelihood = new SpectrumTreeLikelihood(
						spectrumModel, treeModel, siteModel, null, false,
						false, true, false, true);
				
				expected = newSpectrumTreeLikelihood.getLogLikelihood();
				spectrumTreeLikelihood.makeDirty();
				likelihood = spectrumTreeLikelihood.getLogLikelihood();
				System.out.println(likelihood);
				System.out.println(expected);
				
				for (int j = 0; j < taxaCount; j++) {
					System.out.println(treeModel.getTaxon(j));
				}

				SpectrumAlignmentModel defaultSpecturmModer = new SpectrumAlignmentModel(spectrumLength, taxaCount);
				SpectrumAlignmentUtils.compareTwoSpectrumModel(spectrumModel, defaultSpecturmModer);
				
				System.out.println(newSpectrumTreeLikelihood.diagnostic());
				System.out.println( 
						SpectrumTreeLikelihood.compareTwoModels(spectrumTreeLikelihood, newSpectrumTreeLikelihood) );
				
				System.out.println(spectrumTreeLikelihood.diagnostic());
				System.exit(1);
			}
			
		}
	}

	@Test
	public void testUpdateSpectrumMCMC() {

		SiteModel siteModel = MCMCSetupHelper.setupSiteModel();

		int taxaCount = 8;
		Taxon[] taxa = new Taxon[taxaCount];
		for (int i = 0; i < taxa.length; i++) {
			taxa[i] = new Taxon("taxa_"+i);
		}

		TreeModel treeModel = createTreeModel(taxa, 0.01);
		
		String[] seqs = new String[] {
				"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"+
				"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"+
				"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"+
				"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"+
				"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"+
				"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"+
				"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"+
				"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"+
				"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"+
				"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"		
		};

			
		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, taxaCount);
		SpectrumTreeLikelihood spectrumTreeLikelihood = new SpectrumTreeLikelihood(spectrumModel, treeModel,
				siteModel, null, false, false, true, false, false);

		SimpleMCMCOperator op = new DeltaExchangeSingleSpectrumOperator(spectrumModel, 0.1, null);
		double likelihood = spectrumTreeLikelihood.getLogLikelihood();
//		for (int i = 0; i < 1e4; i++) {
//			try {
//				spectrumTreeLikelihood.storeModelState();
//				
//				op.doOperation();
//				double newLikelihood = spectrumTreeLikelihood.getLogLikelihood();
//				double rand = newLikelihood-likelihood;
//
//				if(rand>0.5){
//					spectrumTreeLikelihood.acceptModelState();
//					likelihood = newLikelihood;
//				}
//				else{
//					spectrumTreeLikelihood.restoreModelState();
//				}
//
//			} catch (OperatorFailedException e){
//			}
//			
//		}
//		likelihood = spectrumTreeLikelihood.getLogLikelihood();
//		SpectrumAlignmentModel newSpectrumModel = SpectrumAlignmentModel.duplicateSpectrumAlignmentModel(spectrumModel);
//		SpectrumTreeLikelihood newSpectrumTreeLikelihood = new SpectrumTreeLikelihood(newSpectrumModel, treeModel,
//				siteModel, null, false, false, true, false, false);
//		double expected = newSpectrumTreeLikelihood.getLogLikelihood();
//		likelihood = spectrumTreeLikelihood.getLogLikelihood();
//		assertEquals(expected, likelihood, 0);

//		
//		ShortReadsSpectrumLikelihood srpLikelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
//
//		// CompoundLikelihood
//		HashMap<String, Likelihood> compoundlikelihoods = MCMCSetupHelper.setupCompoundLikelihood(
//				popSize, kappa, coalescent, treeLikelihood, srpLikelihood);
//		Likelihood prior = compoundlikelihoods.get(CompoundLikelihoodParser.PRIOR);
//		Likelihood likelihood = compoundlikelihoods.get(CompoundLikelihoodParser.LIKELIHOOD);
//		Likelihood shortReadLikelihood = compoundlikelihoods.get(ShortReadsSpectrumLikelihood.SHORT_READ_LIKELIHOOD);
//		Likelihood posterior = compoundlikelihoods.get(CompoundLikelihoodParser.POSTERIOR);
//
		// Operators
		OperatorSchedule schedule = new SimpleOperatorSchedule();
		MCMCSetupHelperSpectrum.defalutSpectrumOperators(schedule, spectrumModel);//, freqs, kappa);
//		MCMCSetupHelperSpectrum.defalutTreeOperators(schedule, treeModel);

		// MCLogger
//		MCLogger[] loggers = new MCLogger[1];
//		// log tracer
//		loggers[0] = new MCLogger(new TabDelimitedFormatter(System.out), logInterval, true, logInterval*2);
//		MCMCSetupHelper.addToLogger(loggers[0],
//				posterior, prior, likelihood, shortReadLikelihood,
//				kappa
//				);
		
		// MCMC
		MCMCOptions options = MCMCSetupHelper.setMCMCOptions(10000,1);
		
		MCMC mcmc = new MCMC("mcmc1");
//		mcmc.setShowOperatorAnalysis(true);
//		mcmc.setOperatorAnalysisFile(new File(operatorAnalysisFile));
		
		mcmc.init(options, spectrumTreeLikelihood, schedule, null);
		mcmc.run();
		
		likelihood = spectrumTreeLikelihood.getLogLikelihood();
		SpectrumAlignmentModel newSpectrumModel = SpectrumAlignmentModel.duplicateSpectrumAlignmentModel(spectrumModel);
		SpectrumTreeLikelihood newSpectrumTreeLikelihood = new SpectrumTreeLikelihood(newSpectrumModel, treeModel,
				siteModel, null, false, false, true, false, false);
		double expected = newSpectrumTreeLikelihood.getLogLikelihood();
		likelihood = spectrumTreeLikelihood.getLogLikelihood();
//		assertEquals(expected, likelihood, 0);
		

		// clock model
		Parameter rateParameter = new Parameter.Default(StrictClockBranchRates.RATE, 1e-5, 0, 1);
		StrictClockBranchRates branchRateModel = new StrictClockBranchRates(rateParameter);
	
		Parameter freqs = new Parameter.Default("frequency", new double[]{0.25, 0.25, 0.25, 0.25});
		Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100.0);

		// Sub model
		FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
		HKY hky = new HKY(kappa, f);

		// siteModel
		GammaSiteModel AsiteModel = new GammaSiteModel(hky);
		Parameter mu = new Parameter.Default(
				GammaSiteModelParser.MUTATION_RATE, 1, 0, Double.POSITIVE_INFINITY);
		AsiteModel.setMutationRateParameter(mu);
		
		SpectrumTreeLikelihood anewSpectrumTreeLikelihood = new SpectrumTreeLikelihood(newSpectrumModel, treeModel,
				AsiteModel, null, false, false, true, false, false);
		double aexpected = newSpectrumTreeLikelihood.getLogLikelihood();
		likelihood = anewSpectrumTreeLikelihood.getLogLikelihood();
		assertEquals(aexpected, likelihood, 0);
		
		
	}

	@Test
	public void testUpdateSpectrumAccept() {

		SiteModel siteModel = MCMCSetupHelper.setupSiteModel();

		int taxaCount = 8;
		Taxon[] taxa = new Taxon[taxaCount];
		for (int i = 0; i < taxa.length; i++) {
			taxa[i] = new Taxon("taxa_"+i);
		}

		TreeModel treeModel = createTreeModel(taxa, 0.01);
		
		String[] seqs = new String[] {"ACGTTGCA"};
			
		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, taxaCount);
		SpectrumTreeLikelihood spectrumTreeLikelihood = new SpectrumTreeLikelihood(spectrumModel, treeModel,
				siteModel, null, false, false, true, false, false);

		MCMCOperator op = new DeltaExchangeSingleSpectrumOperator(spectrumModel, 0.1, null);
		double likelihood = spectrumTreeLikelihood.getLogLikelihood();
		
		//Accept only
		for (int i = 0; i < 1e3; i++) {
			try {
				spectrumTreeLikelihood.storeModelState();
				op.operate();
				SpectrumAlignmentModel newSpectrumModel = SpectrumAlignmentModel.duplicateSpectrumAlignmentModel(spectrumModel);
				SpectrumTreeLikelihood newSpectrumTreeLikelihood = new SpectrumTreeLikelihood(newSpectrumModel, treeModel,
						siteModel, null, false, false, true, false, false);
				double expected = newSpectrumTreeLikelihood.getLogLikelihood();
				likelihood = spectrumTreeLikelihood.getLogLikelihood();
				assertEquals(expected, likelihood, 0);
				op.accept(0);
				spectrumTreeLikelihood.acceptModelState();
			}
			catch(OperatorFailedException e){
				op.reject();
			}
		}
	}

	@Test
	public void testUpdateSpectrumRestore() throws OperatorFailedException {

		SiteModel siteModel = MCMCSetupHelper.setupSiteModel();

		int taxaCount = 8;
		Taxon[] taxa = new Taxon[taxaCount];
		for (int i = 0; i < taxa.length; i++) {
			taxa[i] = new Taxon("taxa_"+i);
		}

		TreeModel treeModel = createTreeModel(taxa, 0.01);
		
		String[] seqs = new String[] {"ACGTTGCA"};
		
		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, taxaCount);
		SpectrumTreeLikelihood spectrumTreeLikelihood = new SpectrumTreeLikelihood(spectrumModel, treeModel,
				siteModel, null, false, false, true, false, false);

		SimpleMCMCOperator op = new DeltaExchangeSingleSpectrumOperator(spectrumModel, 0.1, null);
		double likelihood = spectrumTreeLikelihood.getLogLikelihood();
		
		for (int i = 0; i < 1e3; i++) {
			try {
				spectrumTreeLikelihood.storeModelState();
				op.doOperation();
				SpectrumAlignmentModel newSpectrumModel = SpectrumAlignmentModel.duplicateSpectrumAlignmentModel(spectrumModel);
				SpectrumTreeLikelihood newSpectrumTreeLikelihood = new SpectrumTreeLikelihood(newSpectrumModel, treeModel,
						siteModel, null, false, false, true, false, false);
				double expected = newSpectrumTreeLikelihood.getLogLikelihood();
				likelihood = spectrumTreeLikelihood.getLogLikelihood();
				assertEquals(expected, likelihood, 0);
				spectrumTreeLikelihood.restoreModelState();
			}
			catch(OperatorFailedException e){
			}
		}
	}
	
	@Test
	public void testCompareToTreeLikelihood() throws Exception {

		// Sub model
		Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25,
				0.25, 0.25 });
		Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100);

		FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
		HKY hky = new HKY(kappa, f);

		// siteModel
		GammaSiteModel siteModel = new GammaSiteModel(hky);
		Parameter mu = new Parameter.Default(
				GammaSiteModelParser.MUTATION_RATE, 1.0, 0,
				Double.POSITIVE_INFINITY);
		siteModel.setMutationRateParameter(mu);

		// treeLikelihood
		Taxon[] taxa = new Taxon[8];
		for (int i = 0; i < taxa.length; i++) {
			taxa[i] = new Taxon("taxa_"+i);
		}

		TreeModel treeModel = createTreeModel(taxa, 0.01);
		char[] DANChars = new char[]{'A', 'C', 'G', 'T'};
		for (int i = 0; i < 1000; i++) {
			int length = MathUtils.nextInt(100)+10;
			String[] seqs = new String[taxa.length];
			for (int t = 0; t < seqs.length; t++) {
				String s = "";
				for (int c = 0; c < length; c++) {
					s+= DANChars[MathUtils.nextInt(4)];
				}
				seqs[t] = s;
			}
			SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
			AlignmentMapping aMap = new AlignmentMapping(alignment);
			
	//		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 4, 2);
			SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(alignment);
			SpectrumTreeLikelihood spectrumTreeLikelihood = new SpectrumTreeLikelihood(spectrumModel, treeModel,
					siteModel, null, false, false, true, false, false);
	
			HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, alignment);
			TreeLikelihoodExt treeLikelihood = new TreeLikelihoodExt(haplotypeModel, treeModel, 
					siteModel, null, null, false, false, true, false, false);
			
			assertEquals(spectrumTreeLikelihood.getLogLikelihood(), 
					treeLikelihood.getLogLikelihood(), THRESHOLD);
		}

	}
	private static TreeModel createTreeModel(Taxon[] taxa, double height) {
		
		int length = taxa.length;
		if(length == 2){
			return createTreeModel2Taxa(taxa, height);
		}
		else if(length == 4){
			return createTreeModel4Taxa(taxa, height);
		}
		else if(length == 8){
			return createTreeModel8Taxa(taxa, height);
		}
		return null;
	}
	
	private static TreeModel createTreeModel2Taxa(Taxon[] taxa, double height) {

		SimpleNode[] nodes = new SimpleNode[taxa.length*2-1];
		for (int n = 0; n < nodes.length; n++) {
			nodes[n] = new SimpleNode();
		}

		for (int i = 0; i < taxa.length; i++) {
			nodes[i].setTaxon(taxa[i]);
		}
		
		nodes[2].setHeight(height);
		nodes[2].addChild(nodes[0]);
		nodes[2].addChild(nodes[1]);

		SimpleNode root = new SimpleNode();
		root = nodes[2];

		Tree tree = new SimpleTree(root);
		tree.setUnits(Units.Type.YEARS);

		return new TreeModel(tree); // treeModel
	}

	private static TreeModel createTreeModel4Taxa(Taxon[] taxa, double height) {

		SimpleNode[] nodes = new SimpleNode[taxa.length*2-1];
		for (int n = 0; n < nodes.length; n++) {
			nodes[n] = new SimpleNode();
		}

		for (int i = 0; i < taxa.length; i++) {
			nodes[i].setTaxon(taxa[i]);
		}
		
		nodes[4].setHeight(height);
		nodes[4].addChild(nodes[0]);
		nodes[4].addChild(nodes[1]);

		nodes[5].setHeight(height);
		nodes[5].addChild(nodes[2]);
		nodes[5].addChild(nodes[3]);

		nodes[6].setHeight(height*2);
		nodes[6].addChild(nodes[4]);
		nodes[6].addChild(nodes[5]);

		SimpleNode root = new SimpleNode();
		root = nodes[6];
		Tree tree = new SimpleTree(root);
		tree.setUnits(Units.Type.YEARS);

		return new TreeModel(tree); // treeModel
	}

	
	private static TreeModel createTreeModel8Taxa(Taxon[] taxa, double height) {

		SimpleNode[] nodes = new SimpleNode[taxa.length*2-1];
		for (int n = 0; n < nodes.length; n++) {
			nodes[n] = new SimpleNode();
		}

		for (int i = 0; i < taxa.length; i++) {
			nodes[i].setTaxon(taxa[i]);
		}
		
		nodes[8].setHeight(height);
		nodes[8].addChild(nodes[0]);
		nodes[8].addChild(nodes[1]);

		nodes[9].setHeight(height*2);
		nodes[9].addChild(nodes[2]);
		nodes[9].addChild(nodes[3]);

		nodes[10].setHeight(height*3);
		nodes[10].addChild(nodes[8]);
		nodes[10].addChild(nodes[9]);

		
		
		nodes[11].setHeight(height*1);
		nodes[11].addChild(nodes[4]);
		nodes[11].addChild(nodes[5]);

		nodes[12].setHeight(height*2);
		nodes[12].addChild(nodes[6]);
		nodes[12].addChild(nodes[11]);

		nodes[13].setHeight(height*3);
		nodes[13].addChild(nodes[7]);
		nodes[13].addChild(nodes[12]);

		nodes[14].setHeight(height*4);
		nodes[14].addChild(nodes[10]);
		nodes[14].addChild(nodes[13]);

		
		
		SimpleNode root = new SimpleNode();
		root = nodes[14];
//		root.setHeight(0.02);
//		root.addChild(nodes[0]);
//		root.addChild(nodes[3]);

		Tree tree = new SimpleTree(root);
		tree.setUnits(Units.Type.YEARS);

		return new TreeModel(tree); // treeModel
	}

	public void testMultiCreation() throws Exception {
	
			// Sub model
			Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25,
					0.25, 0.25 });
			Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100);
	
			FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
			HKY hky = new HKY(kappa, f);
	
			// siteModel
			GammaSiteModel siteModel = new GammaSiteModel(hky);
			Parameter mu = new Parameter.Default(
					GammaSiteModelParser.MUTATION_RATE, 1.0, 0,
					Double.POSITIVE_INFINITY);
			siteModel.setMutationRateParameter(mu);
	
			// treeLikelihood
			Taxon[] taxa = new Taxon[4];
			for (int i = 0; i < taxa.length; i++) {
				taxa[i] = new Taxon("taxa_"+i);
			}
	
			TreeModel treeModel = createTreeModel(taxa, 0.01);
			
			String[] seqs = new String[] {"ACGT"};
				
			
			int spectrumLength = seqs[0].length();
			SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, 4);
			SpectrumTreeLikelihood spectrumTreeLikelihood = new SpectrumTreeLikelihood(spectrumModel, treeModel,
					siteModel, null, false, false, true, false, false);
	
			SimpleMCMCOperator op = new DeltaExchangeSingleSpectrumOperator(spectrumModel, 0.1, null);
	
	
			for (int i = 0; i < 1e7; i++) {
				if (i%10000 == 0){
					double heapSize = Runtime.getRuntime().totalMemory()/1e6; 
					double heapMaxSize = Runtime.getRuntime().maxMemory()/1e6;
					double heapFreeSize = Runtime.getRuntime().freeMemory()/1e6; 
		
					System.out.println(i +"\t"+ heapSize +"\t"+ heapMaxSize +"\t"+ heapFreeSize);
				}
				try {
					SpectrumTreeLikelihood s = new SpectrumTreeLikelihood(spectrumModel, treeModel,
							siteModel, null, false, false, true, false, false);
					
	
				} catch (IllegalArgumentException e) {
					e.printStackTrace();
	//			} catch (OutOfMemoryError e1) {
	//				System.out.println("ite:"+i);
	//				e1.printStackTrace();
				} catch (Throwable e) {
					e.printStackTrace();
				}
			}
		}

}
