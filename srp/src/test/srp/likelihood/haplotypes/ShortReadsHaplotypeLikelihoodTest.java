package test.srp.likelihood.haplotypes;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.evolution.OperationType;
import srp.evolution.shortreads.ShortReadMapping;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.HaplotypeModel;
import srp.likelihood.haplotypes.ShortReadsHaplotypeLikelihood;
import srp.operator.haplotypes.AbstractHaplotypeOperator;
import srp.operator.haplotypes.BaseSingleOperator;
import srp.operator.haplotypes.BasesMultiOperator;
import srp.operator.haplotypes.ColumnOperator;
import srp.operator.haplotypes.HaplotypeRecombinationOperator;
import srp.operator.haplotypes.HaplotypeSwapSectionOperator;
import srp.operator.spectrum.RecombineSectionSpectrumOperator;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.inference.markovchain.MarkovChain;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCCriterion;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.math.MathUtils;

public class ShortReadsHaplotypeLikelihoodTest {


	public static final double ERROR = ShortReadsHaplotypeLikelihood.ERROR_RATE;
	public static final double NOT_ERROR = ShortReadsHaplotypeLikelihood.NOT_ERROR_RATE;
	private static final double THRESHOLD = MarkovChain.EVALUATION_TEST_THRESHOLD;

	private HaplotypeModel haplotypeModel;
	private ShortReadMapping  srpMap;
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {

	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {


		String dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/";
	
		String trueAlignmentFile = "H4_haplotypes.phyml";

		//		String shortReadFile = "H4_srp.fasta";
		String shortReadFile = "SpectrumTest_50srp_200bp.fasta";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		Alignment shortReads = dataImporter.importShortReads(shortReadFile);
		srpMap = new ShortReadMapping(shortReads);
		
		
//		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		haplotypeModel = new HaplotypeModel(4, srpMap.getLength());
	
	}

	@After
	public void tearDown() throws Exception {
	}

	
	@Test
	public void testConstructor(){
		
		ShortReadsHaplotypeLikelihood likelihood = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		double expectedError = 0.0107;
		assertEquals(expectedError, ShortReadsHaplotypeLikelihood.ERROR_RATE, 0);
		assertEquals(1-expectedError, ShortReadsHaplotypeLikelihood.NOT_ERROR_RATE, 0);
		assertEquals(4, ShortReadsHaplotypeLikelihood.STATE_COUNT, 0);
		assertEquals(19, ShortReadsHaplotypeLikelihood.AMBIGUOUS_STATE_COUNT, 0);
	}
	
	@Test
	public void testCalculateLikelihoodIdentical2() {

		String[] seqs = new String[]{"AACCGGTT"};
	
		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
		
		
//		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
//		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
//		ShortReadMapping srpMap = ShortReadMapping
		HaplotypeModel haplotypeModel = new HaplotypeModel(alignment);
		ShortReadsHaplotypeLikelihood srL = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		
		double logLikelihood = srL.getLogLikelihood();
		double expected = -0.086061253223681313806; //dbinom(0,8,E,log=T)
		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);
		
		seqs = new String[]{
				".....AACCGGTT........",
				".....AACCGGTT........",
				".....AACCGGTT........",
				".....AACCGGTT........",

				};
		alignment = AlignmentUtils.createAlignment(seqs);
		srpMap = new ShortReadMapping(alignment);
		seqs = new String[]{"WWWWWAACCGGTTYYYYYYYY"};
		
		alignment = AlignmentUtils.createAlignment(seqs);
		
		haplotypeModel = new HaplotypeModel(alignment);
		srL = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		
		
		logLikelihood = srL.getLogLikelihood();
		expected = -0.086061253223681313806*4; //dbinom(0,8,E,log=T)
		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);
		
	}



	@Test
		public void testCalculateLikelihoodIdentical() {
	
			String[] seqs = new String[]{"AACCGGTT"};
		
			SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
			ShortReadMapping srpMap = new ShortReadMapping(alignment);
			
			HaplotypeModel spectrumModel = new HaplotypeModel(alignment);
			ShortReadsHaplotypeLikelihood likelihood = new ShortReadsHaplotypeLikelihood(
					spectrumModel, srpMap);
	
			double logLikelihood = likelihood.getLogLikelihood();
			double expected = Math.log(1*NOT_ERROR+0*ERROR)*8;
	//		System.out.println((0.25*NOT_ERROR+0.75*ERROR) +"\t"+ Math.log(1*NOT_ERROR+0*ERROR) );
			assertEquals("0 mismatch",expected, logLikelihood, 1e-10);
		}

	@Test
	public void testCalculateLikelihoodAllDiff() {

		String[] seqs = new String[]{"AACCGGTT"};
		Alignment alignment = AlignmentUtils.createAlignment(seqs);
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
		
		
		String[] haps = new String[]{"TTGGCCAA"};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(hapAlignment);
		ShortReadsHaplotypeLikelihood srL = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);

		double logLikelihood = srL.getLogLikelihood();
		double expected = -36.300092300114215504; //dbinom(0,8,E,log=T)
		assertEquals("8 mismatches",expected, logLikelihood, 1e-10);
		
//		
//		srp.removeSequence(0);
//		srp.addSequence(new Sequence(t, ".....TTGGCCAA......."));
//		aMap = new AlignmentMapping(srp);
//		System.out.println(aMap.getFullSrp(0));
//		System.out.println(aMap.getFragmentSrp(0));
//		srL = new ShortReadLikelihood(aMap, hap);
//		logLikelihood = srL.getLogLikelihoodSelect(BINOMIAL);
//		assertEquals("8 mismatches",expected, logLikelihood, 1e-10);

	}
	
	@Test
	public void testCalculateLikelihood2(){

		String[] srps = new String[]{
//				"AAAACCCCGGGGTTTT",
//				"ACGTACGTACGTACGT",
				"....CCCCGGGG....", // 0/8 + 6/8
				"..AAACGT........", // 3/6 + 2/6
				};
		Alignment alignment = AlignmentUtils.createAlignment(srps);
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
		
		
		String[] haps = new String[]{
				"AAAACCCCGGGGTTTT",
				"ACGTACGTACGTACGT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(hapAlignment);
		ShortReadsHaplotypeLikelihood srL = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);

		double logLikelihood = srL.getLogLikelihood();

		
		double expected = -9.1933572982095146386;
//		E = 0.0107; NE = 1-E 
//		#log(sum(dbinom(c(0,6),8,E))) + log(sum(dbinom(c(3,2),6,E)))
//		log(  E^0*NE^8 + E^6*NE^2    ) + log( E^3*NE^3 + E^2*NE^4     )


		assertEquals(expected, logLikelihood, 1e-10);

		
	}
	
	@Test
	public void testCompoundLikelihood() {


		String[] srps = new String[]{
//				"AAAACCCCGGGGTTTT",
//				"ACGTACGTACGTACGT",
				"....CCCCGGGG....", // 0/8 + 6/8
				"..AAACGT........", // 3/6 + 2/6
				};
		Alignment alignment = AlignmentUtils.createAlignment(srps);
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
		
		
		String[] haps = new String[]{
				"AAAACCCCGGGGTTTT",
				"ACGTACGTACGTACGT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(hapAlignment);
		ShortReadsHaplotypeLikelihood srL = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		
		double logLikelihood = srL.getLogLikelihood();

		
		double expected = -9.1933572982095146386;
//		E = 0.0107; NE = 1-E 
//		#log(sum(dbinom(c(0,6),8,E))) + log(sum(dbinom(c(3,2),6,E)))
//		log(  E^0*NE^8 + E^6*NE^2    ) + log( E^3*NE^3 + E^2*NE^4     )

		assertEquals(expected, logLikelihood, 1e-10);

        List<Likelihood> likelihoods = new ArrayList<Likelihood>();        
        likelihoods.add(srL);
        
        Likelihood CSRLikelihood = new CompoundLikelihood(0, likelihoods);
//        CSRLikelihood.setId(ShortReadLikelihood.NAME);

        logLikelihood = CSRLikelihood.getLogLikelihood();
        assertEquals(expected, logLikelihood, 1e-10);
     
	}	


	@Test
	public void testRandomQuickRun() throws Exception {

		ShortReadsHaplotypeLikelihood srpLikelihood = new ShortReadsHaplotypeLikelihood(
				haplotypeModel, srpMap);

		double likelihood = srpLikelihood.getLogLikelihood();
		double initLikelihood = likelihood;

		MCMCCriterion criterion = new MCMCCriterion();
		double hastingsRatio = 0.0;
		double[] logr = { -Double.MAX_VALUE };

		BaseSingleOperator op = new BaseSingleOperator(haplotypeModel);
		for (int i = 0; i < 1e2; i++) {
			op.doOperation();

			srpLikelihood.storeModelState();
			srpLikelihood.makeDirty();

			double newL = srpLikelihood.getLogLikelihood();

			boolean accept = criterion.accept(likelihood, newL, hastingsRatio,
					logr);
			if (accept) {
				likelihood = newL;
				// srpLikelihood.acceptState();
			} else {
				// haplotypeModel.resreject();
				srpLikelihood.restoreModelState();
			}

		}
		double finalLogLikelihood = srpLikelihood.getLogLikelihood();
		assertTrue(finalLogLikelihood >= initLikelihood);

	}


	@Test
	public void testCalculateLikelihoodFixedSpectrum() {
		String[] seqs = new String[]{
				".AA",
				".AC",
				".GT"
				};
		SimpleAlignment shortreadAlignment = AlignmentUtils.createAlignment(seqs);
		ShortReadMapping srpMap = new ShortReadMapping(shortreadAlignment);
		
		
		seqs = new String[]{"AAA"};
		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		
		HaplotypeModel spectrumModel = new HaplotypeModel(alignment);
		ShortReadsHaplotypeLikelihood likelihood = new ShortReadsHaplotypeLikelihood(spectrumModel, srpMap);

		double[] eachLikelihood = likelihood.unittestMethodGetEachLikelihood();
//		System.out.println(Arrays.toString(eachLikelihood));
		double[] expecteds = new double[]{ 
				0+Math.log(1*NOT_ERROR+0*ERROR)*2,
				0+Math.log(1*NOT_ERROR+0*ERROR)*1+Math.log(0*NOT_ERROR+1*ERROR)*1,
				0+Math.log(1*NOT_ERROR+0*ERROR)*0+Math.log(0*NOT_ERROR+1*ERROR)*2
			};
		assertArrayEquals(expecteds, eachLikelihood, 1e-8);
	}

	@Test
	public void testFullvsMaster() throws Exception {
	
		ShortReadsHaplotypeLikelihood likelihood = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		
		for (int i = 0; i < 1e4; i++) {
			int hapCount = MathUtils.nextInt(7)+3;
			int hapLength = srpMap.getLength();
			haplotypeModel = new HaplotypeModel(hapCount, hapLength);

//				likelihood.makeDirty();
			likelihood = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
			double logLikelihoodFull = likelihood.getLogLikelihood();
//				assertEquals(SpectrumOperation.NONE, likelihood.getOperation());
			
//				likelihood.makeDirty();
			double logLikelihoodMaster = likelihood.calculateSrpLikelihoodFullMaster();
//				assertEquals(SpectrumOperation.NONE, likelihood.getOperation());
			assertEquals(logLikelihoodMaster, logLikelihoodFull, THRESHOLD);

		}
	}
	

	private void assertLikelihoodOperator(HaplotypeModel haplotypeModel,
			OperatorSchedule schedule) {
		
		boolean DEBUG = true;

		int ite = (int) 1e4;
		
		ShortReadsHaplotypeLikelihood likelihood = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		double logLikelihoodOperator = 0;
		double logLikelihoodFull;
		double logLikelihoodMaster;

		for (int i = 0; i < ite; i++) {
//			System.out.println("==== ite: "+i);
			likelihood.storeModelState();
			
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
            
            OperationType expectedSpectrumOperation = ((AbstractHaplotypeOperator) mcmcOperator).getOperationType();
			
			if(operatorSucceeded){
				logLikelihoodOperator = likelihood.getLogLikelihood();
				assertEquals(expectedSpectrumOperation, likelihood.getOperation());
				
				HaplotypeModel haplotypeModelFull = HaplotypeModel.duplicateHaplotypeModel(haplotypeModel);
				ShortReadsHaplotypeLikelihood likelihoodFull = new ShortReadsHaplotypeLikelihood(haplotypeModelFull, srpMap);
				logLikelihoodFull = likelihoodFull.getLogLikelihood();
//				logLikelihoodMaster = likelihood.calculateSrpLikelihoodFullMaster();
				assertEquals(OperationType.NONE, likelihoodFull.getOperation());
				assertEquals(logLikelihoodFull, logLikelihoodOperator, THRESHOLD*100);
				double rand = MathUtils.nextDouble();
				accept = rand>0.5;
			}
			if(accept){
				mcmcOperator.accept(0);
				likelihood.acceptModelState();
			}
			else{
				mcmcOperator.reject();
				likelihood.restoreModelState();
			}
//			if(DEBUG){
//				System.out.println(i);
//			}

		}
		logLikelihoodOperator = likelihood.getLogLikelihood();
		
		HaplotypeModel haplotypeModelFull = HaplotypeModel.duplicateHaplotypeModel(haplotypeModel);
		ShortReadsHaplotypeLikelihood likelihoodFull = new ShortReadsHaplotypeLikelihood(haplotypeModelFull, srpMap);
		logLikelihoodFull = likelihoodFull.getLogLikelihood();
		assertEquals(logLikelihoodFull, logLikelihoodOperator, THRESHOLD);
	}

	@Test
	public void testFullvsSingle() throws Exception {

//		DeltaExchangeSingleSpectrumOperator op = new DeltaExchangeSingleSpectrumOperator(
//				spectrumModel, 0.25, null);
		OperatorSchedule schedule = new SimpleOperatorSchedule();
		MCMCOperator op;
		
//		op = new DirichletAlphaSpectrumOperator(haplotypeModel, 100, null);
//		schedule.addOperator(op);
//		
//		op = new DeltaExchangeSingleSpectrumOperator(haplotypeModel, 0.05, null);
//		schedule.addOperator(op);
		
		op = new BaseSingleOperator(haplotypeModel);
		schedule.addOperator(op);
		
		assertLikelihoodOperator(haplotypeModel, schedule);
	}
	@Test
	public void testFullvsRecombination() throws Exception {

		OperatorSchedule schedule = new SimpleOperatorSchedule();
		MCMCOperator op;
		
		op = new HaplotypeRecombinationOperator(haplotypeModel, 0);
		schedule.addOperator(op);
		op.setWeight(1);
		
		op = new HaplotypeSwapSectionOperator(haplotypeModel, 10, null);
		schedule.addOperator(op);
		op.setWeight(10);
		
		assertLikelihoodOperator(haplotypeModel, schedule);
	}
	
	@Test
	public void testFullvsMulti() throws Exception {

		OperatorSchedule schedule = new SimpleOperatorSchedule();
		MCMCOperator op;
		
		op = new BasesMultiOperator(haplotypeModel, 5, null);
		schedule.addOperator(op);
		
		assertLikelihoodOperator(haplotypeModel, schedule);
	}


	@Test
	public void testFullvsColumnStoreRestore() throws Exception {

		OperatorSchedule schedule = new SimpleOperatorSchedule();
		MCMCOperator op; 

		op = new ColumnOperator(haplotypeModel, haplotypeModel.getHaplotypeCount(), null);
//		op = new DeltaExchangeColumnSpectrumOperator(
//				haplotypeModel, 0.1, null);
		schedule.addOperator(op);
		
		assertLikelihoodOperator(haplotypeModel, schedule);

	}

}
