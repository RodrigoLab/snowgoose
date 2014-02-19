package test.srp.spectrum.likelihood;


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
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import srp.spectrum.likelihood.ShortReadsSpectrumLikelihood;
import srp.spectrum.operator.AbstractSpectrumOperator;
import srp.spectrum.operator.DeltaExchangeColumnSpectrumOperator;
import srp.spectrum.operator.DeltaExchangeMultiSpectrumOperator;
import srp.spectrum.operator.DeltaExchangeSingleSpectrumOperator;
import srp.spectrum.operator.RecombinationSpectrumOperator;
import srp.spectrum.operator.RecombineSectionSpectrumOperator;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.inference.markovchain.MarkovChain;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

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
	private AlignmentMapping aMap;
	@Before
	public void setUp() throws Exception {

//		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "HaplotypeModelTest_10_srp.fasta");
//		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "H4_srp.fasta");
		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "SpectrumTest_50srp_200bp.fasta");
		aMap = new AlignmentMapping(alignment);
		spectrumModel = new SpectrumAlignmentModel(aMap, 4, 2);
	}

	@After
	public void tearDown() throws Exception {
	}


	@Test
	public void testCalculateLikelihoodIdentical() {

		String[] seqs = new String[]{"AACCGGTT"};
	
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
//		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		
		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, alignment);
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
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
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		seqs = new String[]{"AAA"};
		SimpleAlignment alignment = AlignmentUtils.createAlignment(seqs);
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, alignment);
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
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
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
			
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 1, 0);
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
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
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 1, 2);
		Spectrum spectrum = spectrumModel.getSpectrum(0);
		for (int i = 0; i < spectrum.getLength(); i++) {
			double[] freqs = new double[]{1-(0.1*i*3), 0.1*i, 0.1*i, 0.1*i};
			spectrum.resetFrequencies(i, freqs);
//			System.out.println("SITE: "+i +"\t"+  Arrays.toString(spectrum.getFrequencies(i)));
		}
//		spectrumModel.setSpectrum(0, spectrum);
		
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
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
		

	}
	@Test
	public void testFullvsMaster() throws Exception {
	
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		for (int i = 0; i < 1e4; i++) {
			int noSpectrum = MathUtils.nextInt(7)+3;
			spectrumModel = new SpectrumAlignmentModel(aMap, noSpectrum, 2);
//				likelihood.makeDirty();
			likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
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

	private static void assertLikelihoodOperator(SpectrumAlignmentModel spectrumModel,
			AbstractSpectrumOperator op) {
		
		int ite = (int) 1e4;
		SpectrumOperation expectedSpectrumOperation = op.getSpectrumOperation();
		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		double logLikelihoodOperator;
		double logLikelihoodFull;

		for (int i = 0; i < ite; i++) {
			try {
				likelihood.storeModelState();
				
				op.doOperation();

				logLikelihoodOperator = likelihood.getLogLikelihood();
				assertEquals(expectedSpectrumOperation, likelihood.getOperation());
				
				SpectrumAlignmentModel spectrumModelFull = SpectrumAlignmentModel.duplicateSpectrumAlignmentModel(spectrumModel);
				ShortReadsSpectrumLikelihood likelihoodFull = new ShortReadsSpectrumLikelihood(spectrumModelFull);
				logLikelihoodFull = likelihoodFull.getLogLikelihood();
				assertEquals(SpectrumOperation.NONE, likelihoodFull.getOperation());
				assertEquals(logLikelihoodFull, logLikelihoodOperator, THRESHOLD); 

	
				double rand = MathUtils.nextDouble();
				if(rand>0.5){
					likelihood.acceptModelState();
				}
				else{
					likelihood.restoreModelState();
				}

			} catch (OperatorFailedException e) {
			}
		}
	}

	@Test
	public void testFullvsSingleStoreRestore() throws Exception {

		DeltaExchangeSingleSpectrumOperator op = new DeltaExchangeSingleSpectrumOperator(
				spectrumModel, 0.25, null);
		assertLikelihoodOperator(spectrumModel, op);
	}

	@Test
	public void testFullvsMultiStoreRestore() throws Exception {

		DeltaExchangeMultiSpectrumOperator op = new DeltaExchangeMultiSpectrumOperator(
				spectrumModel, 0.1, 10, null);
		assertLikelihoodOperator(spectrumModel, op);
	}

	@Test
	public void testFullvsColumnStoreRestore() throws Exception {

		DeltaExchangeColumnSpectrumOperator op = new DeltaExchangeColumnSpectrumOperator(
				spectrumModel, 0.1, null);
		assertLikelihoodOperator(spectrumModel, op);

	}

	@Test
	public void testFullvsRecombinationStoreRestore() throws Exception {

		RecombinationSpectrumOperator op = new RecombinationSpectrumOperator(spectrumModel);
		assertLikelihoodOperator(spectrumModel, op);
	}

	@Test
	public void testFullvsSwapSectionStoreRestore() throws Exception {

		RecombineSectionSpectrumOperator op = new RecombineSectionSpectrumOperator(
				spectrumModel, 10, null);
		assertLikelihoodOperator(spectrumModel, op);
	}

	@Test
	public void testCalculateStateLogLikelihoodFlat() throws Exception{
		
		ShortReadsSpectrumLikelihood srpLikelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		Method method = ShortReadsSpectrumLikelihood.class.getDeclaredMethod("caluclateStateLogLikelihood", double.class);
		method.setAccessible(true);
		
		double[] frequencies = new double[21];
		for (int i = 0; i < frequencies.length; i++) {
			frequencies[i]= i/20.0;
		}
		frequencies[0] = 0.001;
		frequencies[20] = 0.999;
		
		srpLikelihood.setDistType("flat");
		double[] expecteds = new double[] { -4.4499971717799,
				-2.81959647584712, -2.22045226345552, -1.84839331479565,
				-1.57784235084435, -1.36512012588848, -1.18980694885304,
				-1.04069249808234, -0.910954992284636, -0.796132740880369,
				-0.693147180559945, -0.599784350155737, -0.514398666152773,
				-0.43573361212845, -0.362807998442184, -0.294840969650261,
				-0.231200924941943, -0.17136974738986, -0.114917146243339,
				-0.0614818641436684, -0.0117473304903028 };
		for (int i = 0; i < frequencies.length; i++) {
			double log = (Double) method.invoke(srpLikelihood, frequencies[i]);
			assertEquals("fail at "+i, expecteds[i], log, THRESHOLD);
		}
		
	}
	@Test
	public void testCalculateStateLogLikelihoodBetaMean() throws Exception{
		
		ShortReadsSpectrumLikelihood srpLikelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		Method method = ShortReadsSpectrumLikelihood.class.getDeclaredMethod("caluclateStateLogLikelihood", double.class);
		method.setAccessible(true);
		
		double[] frequencies = new double[21];
		for (int i = 0; i < frequencies.length; i++) {
			frequencies[i]= i/20.0;
		}
		frequencies[0] = 0.001;
		frequencies[20] = 0.999;
	
		srpLikelihood.setDistType("betaMean");
		double[] expecteds = new double[] { -4.4627970966445,
				-4.45490108164499, -4.4088290544743, -4.35662071231919,
				-4.29972297213143, -4.23826264916909, -4.17195884206474,
				-4.10029323748722, -4.02253577268687, -3.93771569591255,
				-3.84455269254967, -3.74133935333329, -3.62574761789222,
				-3.49450146814089, -3.34279315588983, -3.16316066347848,
				-2.9430953103395, -2.65914011951639, -2.25862508309307,
				-1.57347309663271, 2.29615312974098 };
		for (int i = 0; i < 21; i++) {
			double log = (Double) method.invoke(srpLikelihood, frequencies[i]);
			assertEquals("fail at "+i, expecteds[i], log, THRESHOLD);
		}
	}
	
	@Test
	public void testCalculateStateLogLikelihoodBetaMode() throws Exception{
		
		ShortReadsSpectrumLikelihood srpLikelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		Method method = ShortReadsSpectrumLikelihood.class.getDeclaredMethod("caluclateStateLogLikelihood", double.class);
		method.setAccessible(true);
		
		double[] frequencies = new double[21];
		for (int i = 0; i < frequencies.length; i++) {
			frequencies[i]= i/20.0;
		}
		frequencies[0] = 0.001;
		frequencies[20] = 0.999;
	
		srpLikelihood.setDistType("betaMode");
		double[] expecteds = new double[] { -6.13013650123304,
				-2.26051027485935, -1.57535828839899, -1.17484325197567,
				-0.890888061152559, -0.670822708013584, -0.491190215602235,
				-0.339481903351173, -0.20823575359984, -0.0926440181587687,
				0.0105693210576128, 0.103732324420493, 0.18855240119481,
				0.266309865995157, 0.337975470572683, 0.404279277677028,
				0.465739600639366, 0.522637340827131, 0.574845682982234,
				0.620917710152927, 0.628813725152435 };
		for (int i = 0; i < 21; i++) {
			double log = (Double) method.invoke(srpLikelihood, frequencies[i]);
			assertEquals("fail at " + i, expecteds[i], log, THRESHOLD);
		}
	}
	

	@Test
	public void testCalculateStateLogLikelihoodGTest() throws Exception{
		
		ShortReadsSpectrumLikelihood srpLikelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		Method method = ShortReadsSpectrumLikelihood.class.getDeclaredMethod("caluclateStateLogLikelihood", double.class);
		method.setAccessible(true);
		
		double[] frequencies = new double[21];
		for (int i = 0; i < frequencies.length; i++) {
			frequencies[i]= i/20.0;
		}
		frequencies[0] = 0.001;
		frequencies[20] = 0.999;
	
		srpLikelihood.setDistType("gTest");
		double[] expecteds = new double[] { -6.54540700719581,
				-6.08520556397393, -5.68744195480934, -5.31842484136728,
				-4.96806776933106, -4.63151576637113, -4.30581779805705,
				-3.98890554821967, -3.67915590009021, -3.37515321539537,
				-3.07552433951212, -2.77878684442851, -2.48316519890541,
				-2.18631510268386, -1.88484260502715, -1.57335518321388,
				-1.24232964427087, -0.872435902109154, -0.414878166309604,
				0.323417013652523, 1.18179430054454 };
		for (int i = 0; i < 21; i++) {
			double log = (Double) method.invoke(srpLikelihood, frequencies[i]);
			assertEquals("fail at " + i, expecteds[i], log, THRESHOLD);
		}
	}
	

	@Test
	public void testCalculateStateLogLikelihoodChisq() throws Exception{
		
		ShortReadsSpectrumLikelihood srpLikelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		Method method = ShortReadsSpectrumLikelihood.class.getDeclaredMethod("caluclateStateLogLikelihood", double.class);
		method.setAccessible(true);
		
		double[] frequencies = new double[21];
		for (int i = 0; i < frequencies.length; i++) {
			frequencies[i]= i/20.0;
		}
		frequencies[0] = 0.001;
		frequencies[20] = 0.999;
	
		srpLikelihood.setDistType("chisq");
		double[] expecteds = new double[] { -49.3168653997597,
				-44.8046134842287, -40.4312738560728, -36.2909399823289,
				-32.3832226986611, -28.7076563643503, -25.2636774106532,
				-22.0505947966424, -19.0675484079555, -16.3134489648306,
				-13.786888606292, -11.4860031036221, -9.40825045103395,
				-7.55003641111874, -5.90603939335547, -4.46788851746636,
				-3.22127082107978, -2.13850678309654, -1.15398942877333,
				-0.0294953974800052, 1.43811198038661 };
		for (int i = 0; i < 21; i++) {
			double log = (Double) method.invoke(srpLikelihood, frequencies[i]);
			assertEquals("fail at " + i, expecteds[i], log, THRESHOLD);
		}
	}
	
}
