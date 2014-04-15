package test.srp.likelihood.haplotypes;


import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.evolution.OperationType;
import srp.evolution.shortreads.ShortReadMapping;
import srp.evolution.spectrum.Spectrum;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.HaplotypeModel;
import srp.likelihood.haplotypes.ShortReadsHaplotypeLikelihood;
import srp.operator.haplotypes.AbstractHaplotypeOperator;
import srp.operator.haplotypes.BaseSingleOperator;
import srp.operator.haplotypes.BasesMultiOperator;
import srp.operator.spectrum.DeltaExchangeColumnSpectrumOperator;
import srp.operator.spectrum.DeltaExchangeMultiSpectrumOperator;
import srp.operator.spectrum.DirichletSpectrumOperator;
import srp.operator.spectrum.RecombinationSpectrumOperator;
import srp.operator.spectrum.RecombineSectionSpectrumOperator;
import srp.operator.spectrum.SwapMultiSpectrumOperator;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.inference.markovchain.MarkovChain;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.math.MathUtils;

/*
 * reflection example
 * 	
ShortReadsHaplotypeLikelihood srpLikelihood = new ShortReadsHaplotypeLikelihood(spectrumModel);
Method method = ShortReadsHaplotypeLikelihood.class.getDeclaredMethod("caluclateStateLogLikelihood", double.class);
method.setAccessible(true);
double log = (Double) method.invoke(srpLikelihood, 0.5);

 */

public class ShortReadsHaplotypeLikelihoodOperatorTest {

	public static final double ERROR = ShortReadsHaplotypeLikelihood.ERROR_RATE;
	public static final double NOT_ERROR = ShortReadsHaplotypeLikelihood.NOT_ERROR_RATE;
	private static final double THRESHOLD = MarkovChain.EVALUATION_TEST_THRESHOLD;
//	private static final double EVALUATION_TEST_THRESHOLD = 1e-8;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	private HaplotypeModel haplotypeModel;
	private ShortReadMapping srpMap;

	@Before
	public void setUp() throws Exception {

//		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "HaplotypeModelTest_10_srp.fasta");
//		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "H4_srp.fasta");
		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "SpectrumTest_50srp_200bp.fasta");
		srpMap = new ShortReadMapping(alignment);
		haplotypeModel = new HaplotypeModel(5, srpMap.getLength());
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

		int ite = (int) 1e2;
		
		ShortReadsHaplotypeLikelihood likelihood = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		double logLikelihoodOperator = 0;
		double logLikelihoodFull;
		double logLikelihoodMaster;

		for (int i = 0; i < ite; i++) {
			System.out.println("================== ite: "+i);
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
				assertEquals(logLikelihoodFull, logLikelihoodOperator, THRESHOLD);
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
//		logLikelihoodOperator = likelihood.getLogLikelihood();
//		
//		HaplotypeModel haplotypeModelFull = HaplotypeModel.duplicateHaplotypeModel(haplotypeModel);
//		ShortReadsHaplotypeLikelihood likelihoodFull = new ShortReadsHaplotypeLikelihood(haplotypeModelFull, srpMap);
//		logLikelihoodFull = likelihoodFull.getLogLikelihood();
//		assertEquals(logLikelihoodFull, logLikelihoodOperator, THRESHOLD);
	}

	@Test
	public void testFullvsSingleStoreRestore() throws Exception {

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
	public void testFullvsMultiStoreRestore() throws Exception {

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

		op = new DeltaExchangeColumnSpectrumOperator(
				haplotypeModel, 0.1, null);
		schedule.addOperator(op);
		
		assertLikelihoodOperator(haplotypeModel, schedule);

	}

	@Test
	public void testFullvsRecombinationStoreRestore() throws Exception {

		OperatorSchedule schedule = new SimpleOperatorSchedule();
		MCMCOperator op;
		
		op = new RecombinationSpectrumOperator(haplotypeModel);
		schedule.addOperator(op);
		
		op = new RecombineSectionSpectrumOperator(haplotypeModel, 2, null);
		schedule.addOperator(op);
		
		assertLikelihoodOperator(haplotypeModel, schedule);
	}



	
	
}
