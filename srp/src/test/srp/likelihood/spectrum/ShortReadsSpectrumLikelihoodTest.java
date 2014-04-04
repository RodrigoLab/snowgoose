package test.srp.likelihood.spectrum;


import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.lang.reflect.Method;
import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.evolution.OperationType;
import srp.evolution.shortreads.AlignmentMapping;
import srp.evolution.shortreads.ShortReadMapping;
import srp.evolution.spectrum.Spectrum;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import srp.evolution.spectrum.SpectraParameter.SpectraType;
import srp.haplotypes.AlignmentUtils;
import srp.likelihood.spectrum.ShortReadsSpectrumLikelihood;
import srp.operator.spectrum.AbstractSpectrumOperator;
import srp.operator.spectrum.DeltaExchangeColumnSpectrumOperator;
import srp.operator.spectrum.DeltaExchangeMultiSpectrumOperator;
import srp.operator.spectrum.DeltaExchangeSingleSpectrumOperator;
import srp.operator.spectrum.DirichletAlphaSpectrumOperator;
import srp.operator.spectrum.DirichletSpectrumOperator;
import srp.operator.spectrum.RecombinationSpectrumOperator;
import srp.operator.spectrum.RecombineSectionSpectrumOperator;
import srp.operator.spectrum.SwapMultiSpectrumOperator;
import srp.operator.spectrum.SwapSingleSpectrumOperator;
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
ShortReadsSpectrumLikelihood srpLikelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
Method method = ShortReadsSpectrumLikelihood.class.getDeclaredMethod("caluclateStateLogLikelihood", double.class);
method.setAccessible(true);
double log = (Double) method.invoke(srpLikelihood, 0.5);

 */

public class ShortReadsSpectrumLikelihoodTest {

	public static final double ERROR = ShortReadsSpectrumLikelihood.ERROR_RATE;
	public static final double NOT_ERROR = ShortReadsSpectrumLikelihood.NOT_ERROR_RATE;
	private static final double THRESHOLD = MarkovChain.EVALUATION_TEST_THRESHOLD;
//	private static final double EVALUATION_TEST_THRESHOLD = 1e-8;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	private SpectrumAlignmentModel spectrumModel;
	private ShortReadMapping srpMap;

	@Before
	public void setUp() throws Exception {

//		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "HaplotypeModelTest_10_srp.fasta");
//		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "H4_srp.fasta");
		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "SpectrumTest_50srp_200bp.fasta");
		srpMap = new ShortReadMapping(alignment);
		spectrumModel = new SpectrumAlignmentModel(alignment.getSiteCount(),  4);
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void testConstructor(){
		
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap);
		double expectedError = 0.0107;
		assertEquals(expectedError, ShortReadsSpectrumLikelihood.ERROR_RATE, 0);
		assertEquals(1-expectedError, ShortReadsSpectrumLikelihood.NOT_ERROR_RATE, 0);
		assertEquals(4, ShortReadsSpectrumLikelihood.STATE_COUNT, 0);
		assertEquals(19, ShortReadsSpectrumLikelihood.AMBIGUOUS_STATE_COUNT, 0);
	}


	@Test
	public void testCalculateLikelihoodIdentical() {

		String[] seqs = new String[]{"AACCGGTT"};
	
		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(
				alignment);
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(
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
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(alignment);
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap);

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
	public void testCalculateLikelihoodSpectrum() {
		String[] seqs = new String[]{
				".AA",
				".AC",
				".GT"
				};
		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
		
		int spectrumLength = 3;
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, 1, SpectraType.EQUAL);
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap);

		double[] eachLikelihood = likelihood.unittestMethodGetEachLikelihood();
		double[] expecteds = new double[]{ 
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*2,
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*1+Math.log(0.25*NOT_ERROR+0.75*ERROR)*1,
				0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*0+Math.log(0.25*NOT_ERROR+0.75*ERROR)*2
			};
		assertArrayEquals(expecteds, eachLikelihood, 1e-8);
		
//		logLikelihood = likelihood .getLogLikelihood();
//		expected = -0.086061253223681313806*4; //dbinom(0,8,E,log=T)
//		assertEquals("0 mismatch",expected, logLikelihood, 1e-10);

//		double[] expecteds = new double[]{ 
//				0+Math.log((0.25*NOT_ERROR+0.75*ERROR)*2),
//				0+Math.log((0.25*NOT_ERROR+0.75*ERROR)*1)+Math.log((0.75*NOT_ERROR+0.25*ERROR)*1),
//				0+Math.log((0.25*NOT_ERROR+0.75*ERROR)*1)+Math.log((0.75*NOT_ERROR+0.25*ERROR)*1)
//			};

	}



	@Test
	public void testCalculateLikelihoodCustomSpectrum() {
		String[] seqs = new String[]{
				"AAAC",
				"AACT",
				"ACGT"
				};
		
		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
		
		int spectrumLength = seqs[0].length();
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, 1);
		Spectrum spectrum = spectrumModel.getSpectrum(0);
		for (int i = 0; i < spectrum.getLength(); i++) {
			double[] freqs = new double[]{1-(0.1*i*3), 0.1*i, 0.1*i, 0.1*i};
			spectrum.resetSpectra(i, freqs);
//			System.out.println("SITE: "+i +"\t"+  Arrays.toString(spectrum.getFrequencies(i)));
		}
//		spectrumModel.setSpectrum(0, spectrum);
		
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap);

		double[] eachLikelihood = likelihood.unittestMethodGetEachLikelihood();
//		System.out.println(Arrays.toString(eachLikelihood));
		//Site1: 1, 0, 0, 0
		//Site2: 0.7, 0.1, 0.1, 0.1
		//Site3: 0.4, 0.2, 0.2, 0.2
		//Site4: 0.1, 0.3, 0.3, 0.3
		double[] expecteds = new double[] {
				// MMMD
				Math.log((1 * NOT_ERROR + 0 * ERROR)
						* (0.7 * NOT_ERROR + 0.3 * ERROR)
						* (0.4 * NOT_ERROR + 0.6 * ERROR)
						* (0.3 * NOT_ERROR + 0.7 * ERROR)),
				// MMDD
				Math.log((1 * NOT_ERROR + 0 * ERROR)
						* (0.7 * NOT_ERROR + 0.3 * ERROR)
						* (0.2 * NOT_ERROR + 0.8 * ERROR)
						* (0.3 * NOT_ERROR + 0.7 * ERROR)),

				// MDDD
				Math.log((1 * NOT_ERROR + 0 * ERROR)
						* (0.1 * NOT_ERROR + 0.9 * ERROR)
						* (0.2 * NOT_ERROR + 0.8 * ERROR)
						* (0.3 * NOT_ERROR + 0.7 * ERROR)) 
		};
		assertArrayEquals(expecteds, eachLikelihood, 1e-8);
		
		//
		seqs = new String[]{
				"AACC",
				"AACC",
				"GGTT",
				"GGTT"
				};
		
		alignment = AlignmentUtils.createAlignment(seqs);
		srpMap = new ShortReadMapping(alignment);
		
		spectrumLength = seqs[0].length();
		spectrumModel = new SpectrumAlignmentModel(spectrumLength, 2);

		spectrum = spectrumModel.getSpectrum(0);
		double[] freqs = new double[]{0.5, 0, 0.5, 0};
		spectrum.resetSpectra(0, freqs);
		spectrum.resetSpectra(1, freqs);
		freqs = new double[]{0, 0.5, 0, 0.5};
		spectrum.resetSpectra(2, freqs);
		spectrum.resetSpectra(3, freqs);
		
		spectrum = spectrumModel.getSpectrum(1);
		freqs = new double[]{0.5, 0, 0.5, 0};
		spectrum.resetSpectra(0, freqs);
		spectrum.resetSpectra(1, freqs);
		freqs = new double[]{0, 0.5, 0, 0.5};
		spectrum.resetSpectra(2, freqs);
		spectrum.resetSpectra(3, freqs);
		
		String[] trueSeq = new String[]{
				"AACC",
				"GGTT"
		};
		SpectrumAlignmentModel spectrumModel2 = new SpectrumAlignmentModel(spectrumLength, 2);
//		SpectrumAlignmentModel spectrumModel2 = new SpectrumAlignmentModel(AlignmentUtils.createAlignment(trueSeq) );
		spectrum = spectrumModel2.getSpectrum(0);
		freqs = new double[]{1, 0, 0, 0};
		spectrum.resetSpectra(0, freqs);
		spectrum.resetSpectra(1, freqs);
		freqs = new double[]{0, 1, 0, 0};
		spectrum.resetSpectra(2, freqs);
		spectrum.resetSpectra(3, freqs);
		
		spectrum = spectrumModel2.getSpectrum(1);
		freqs = new double[]{0, 0, 1, 0};
		spectrum.resetSpectra(0, freqs);
		spectrum.resetSpectra(1, freqs);
		freqs = new double[]{0, 0, 0, 1};
		spectrum.resetSpectra(2, freqs);
		spectrum.resetSpectra(3, freqs);
//			System.out.println("SITE: "+i +"\t"+  Arrays.toString(spectrum.getFrequencies(i)));
		
//		spectrumModel.setSpectrum(0, spectrum);
		
		likelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap);

		eachLikelihood = likelihood.unittestMethodGetEachLikelihood();
		double logLikelihood = likelihood.getLogLikelihood();
		System.out.println(logLikelihood);
		System.out.println(Arrays.toString(eachLikelihood));

		
		likelihood = new ShortReadsSpectrumLikelihood(spectrumModel2, srpMap);
		eachLikelihood = likelihood.unittestMethodGetEachLikelihood();
		logLikelihood = likelihood.getLogLikelihood();
		System.out.println(logLikelihood);
		
		System.out.println(Arrays.toString(eachLikelihood));
		//Site1: 1, 0, 0, 0
		//Site2: 0.7, 0.1, 0.1, 0.1
		//Site3: 0.4, 0.2, 0.2, 0.2
		//Site4: 0.1, 0.3, 0.3, 0.3
//		double[] expecteds = new double[] {
//				// MMMD
//				Math.log((1 * NOT_ERROR + 0 * ERROR)
//						* (0.7 * NOT_ERROR + 0.3 * ERROR)
//						* (0.4 * NOT_ERROR + 0.6 * ERROR)
//						* (0.3 * NOT_ERROR + 0.7 * ERROR)),
//				// MMDD
//				Math.log((1 * NOT_ERROR + 0 * ERROR)
//						* (0.7 * NOT_ERROR + 0.3 * ERROR)
//						* (0.2 * NOT_ERROR + 0.8 * ERROR)
//						* (0.3 * NOT_ERROR + 0.7 * ERROR)),
//
//				// MDDD
//				Math.log((1 * NOT_ERROR + 0 * ERROR)
//						* (0.1 * NOT_ERROR + 0.9 * ERROR)
//						* (0.2 * NOT_ERROR + 0.8 * ERROR)
//						* (0.3 * NOT_ERROR + 0.7 * ERROR)) 
//		};
//		assertArrayEquals(expecteds, eachLikelihood, 1e-8);
		

	}
	@Test
	public void testFullvsMaster() throws Exception {
	
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap);
		
		for (int i = 0; i < 1e3; i++) {
			int noSpectrum = MathUtils.nextInt(7)+3;
			int spectrumLength = srpMap.getLength();
			spectrumModel = new SpectrumAlignmentModel(spectrumLength, noSpectrum);
//				likelihood.makeDirty();
			likelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap);
			double logLikelihoodFull = likelihood.getLogLikelihood();
//				assertEquals(SpectrumOperation.NONE, likelihood.getOperation());
			
//				likelihood.makeDirty();
			double logLikelihoodMaster = likelihood.calculateSrpLikelihoodFullMaster();
//				assertEquals(SpectrumOperation.NONE, likelihood.getOperation());
			assertEquals(logLikelihoodMaster, logLikelihoodFull, THRESHOLD);

		}
	}
	
//	
//	@Test
//	public void testFullvsSingle() throws Exception {
//	
//		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
//		
////		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
//		DeltaExchangeSingleSpectrumOperator op = new DeltaExchangeSingleSpectrumOperator(spectrumModel, 0.25, null);
//		
//		for (int i = 0; i < 1e4; i++) {
//			try {
//				op.doOperation();
//				double logLikelihoodSingle = likelihood.getLogLikelihood();
//				assertEquals(SpectrumOperation.DELTA_SINGLE, likelihood.getOperation());
//				
//				likelihood.makeDirty();
//				double logLikelihoodFull = likelihood.getLogLikelihood();
//				assertEquals(SpectrumOperation.FULL, likelihood.getOperation());
//				assertEquals(logLikelihoodFull, logLikelihoodSingle, THRESHOLD);
//				
//			} catch (Exception e) {
//			}
//		}
//	}

	private void assertLikelihoodOperator(SpectrumAlignmentModel spectrumModel,
			OperatorSchedule schedule) {
		
		boolean DEBUG = true;

		int ite = (int) 1e4;
		
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap);
		double logLikelihoodOperator;
		double logLikelihoodFull;

		for (int i = 0; i < ite; i++) {
			
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
            
            OperationType expectedSpectrumOperation = ((AbstractSpectrumOperator) mcmcOperator).getSpectrumOperation();
			
			if(operatorSucceeded){
				logLikelihoodOperator = likelihood.getLogLikelihood();
				assertEquals(expectedSpectrumOperation, likelihood.getOperation());
				
				SpectrumAlignmentModel spectrumModelFull = SpectrumAlignmentModel.duplicateSpectrumAlignmentModel(spectrumModel);
				ShortReadsSpectrumLikelihood likelihoodFull = new ShortReadsSpectrumLikelihood(spectrumModelFull, srpMap);
				logLikelihoodFull = likelihoodFull.getLogLikelihood();
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
	}

	@Test
	public void testFullvsSingleStoreRestore() throws Exception {

//		DeltaExchangeSingleSpectrumOperator op = new DeltaExchangeSingleSpectrumOperator(
//				spectrumModel, 0.25, null);
		OperatorSchedule schedule = new SimpleOperatorSchedule();
		MCMCOperator op;
		
		op = new DirichletAlphaSpectrumOperator(spectrumModel, 100, null);
		schedule.addOperator(op);
		
		op = new DeltaExchangeSingleSpectrumOperator(spectrumModel, 0.05, null);
		schedule.addOperator(op);
		
		op = new SwapSingleSpectrumOperator(spectrumModel, true);
		schedule.addOperator(op);
		
		assertLikelihoodOperator(spectrumModel, schedule);
	}

	@Test
	public void testFullvsMultiStoreRestore() throws Exception {

		OperatorSchedule schedule = new SimpleOperatorSchedule();
		MCMCOperator op;
		
		op = new DirichletSpectrumOperator(spectrumModel, 5, 100, null);
		schedule.addOperator(op);
		
		op = new DeltaExchangeMultiSpectrumOperator(spectrumModel, 3, 0.1, null);
		schedule.addOperator(op);
		
		op = new SwapMultiSpectrumOperator(spectrumModel, 3, true, null);
		schedule.addOperator(op);
		
		assertLikelihoodOperator(spectrumModel, schedule);
	}

	@Test
	public void testFullvsColumnStoreRestore() throws Exception {

		OperatorSchedule schedule = new SimpleOperatorSchedule();
		MCMCOperator op; 

		op = new DeltaExchangeColumnSpectrumOperator(
				spectrumModel, 0.1, null);
		schedule.addOperator(op);
		
		assertLikelihoodOperator(spectrumModel, schedule);

	}

	@Test
	public void testFullvsRecombinationStoreRestore() throws Exception {

		OperatorSchedule schedule = new SimpleOperatorSchedule();
		MCMCOperator op;
		
		op = new RecombinationSpectrumOperator(spectrumModel);
		schedule.addOperator(op);
		
		op = new RecombineSectionSpectrumOperator(spectrumModel, 2, null);
		schedule.addOperator(op);
		
		assertLikelihoodOperator(spectrumModel, schedule);
	}



	
	
}
