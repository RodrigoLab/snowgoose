package test.srp.spectrum.likelihood.stateLikelihood;

import static org.junit.Assert.*;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.haplotypes.AlignmentMapping;
import srp.spectrum.SpectraParameter;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectraParameter.SpectraType;
import srp.spectrum.likelihood.ShortReadsSpectrumLikelihood;
import srp.spectrum.likelihood.stateLikelihood.BetaMeanStateLikelihood;
import srp.spectrum.likelihood.stateLikelihood.BetaModeStateLikelihood;
import srp.spectrum.likelihood.stateLikelihood.ChisqStateLikelihood;
import srp.spectrum.likelihood.stateLikelihood.GTestStateLikelihood;
import srp.spectrum.likelihood.stateLikelihood.ProbabilityStateLikelihood;
import srp.spectrum.likelihood.stateLikelihood.StateLikelihood;
import dr.evolution.alignment.Alignment;
import dr.inference.markovchain.MarkovChain;

public class StateLikelihoodTest {

	private static final double THRESHOLD = MarkovChain.EVALUATION_TEST_THRESHOLD;
//	private SpectrumAlignmentModel spectrumModel;
//	private AlignmentMapping aMap;
	double[] frequencies = new double[21];
	

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}



	@Before
	public void setUp() throws Exception {

		for (int i = 0; i < frequencies.length; i++) {
			frequencies[i]= i/20.0;
		}
		frequencies[0] = 0.001;
		frequencies[20] = 0.999;
	}
	
	@After
	public void tearDown() throws Exception {
	}

	
	@Test
	public void testCalculateStatesLogLikelihood() throws Exception {
		
		double spectraFrequencies[][] = new double[5][4];
		for (int i = 0; i < spectraFrequencies.length; i++) {
			System.arraycopy(frequencies, 4*i, spectraFrequencies[i], 0, 4);
			System.out.println(Arrays.toString(spectraFrequencies[i]));
			
			
		}
		StateLikelihood stateLikelihood = new ProbabilityStateLikelihood();
		double[] statesLogLikelihood = new double[4];
		double[] expectedLogLikelihood = new double[4];
		double[] expecteds = new double[] { -4.4499971717799,
				-2.81959647584712, -2.22045226345552, -1.84839331479565,
				-1.57784235084435, -1.36512012588848, -1.18980694885304,
				-1.04069249808234, -0.910954992284636, -0.796132740880369,
				-0.693147180559945, -0.599784350155737, -0.514398666152773,
				-0.43573361212845, -0.362807998442184, -0.294840969650261,
				-0.231200924941943, -0.17136974738986, -0.114917146243339,
				-0.0614818641436684};
		double[] expectedLast = new double[]{-0.0117473304903028, -0.0117473304903028, -0.0117473304903028, -0.0117473304903028}; 
		for (int i = 0; i < spectraFrequencies.length; i++) {
			SpectraParameter sp = new SpectraParameter(SpectraType.EQUAL);
			for (int j = 0; j < spectraFrequencies[i].length; j++) {
				sp.setFrequency(j, spectraFrequencies[i][j]);
			}
			
			stateLikelihood.calculateStatesLogLikelihood(sp, statesLogLikelihood);
			System.arraycopy(expecteds, 4*i, expectedLogLikelihood, 0, 4);
			assertArrayEquals(expectedLogLikelihood, statesLogLikelihood, THRESHOLD);
			
			sp.storeState();
			for (int j = 0; j < spectraFrequencies[i].length; j++) {
				sp.setFrequency(j, 0.999);
			}
			stateLikelihood.calculateStoredStatesLogLikelihood(sp, statesLogLikelihood);
			assertArrayEquals(expectedLogLikelihood, statesLogLikelihood, THRESHOLD);
			
			stateLikelihood.calculateStatesLogLikelihood(sp, statesLogLikelihood);
			assertArrayEquals(expectedLast, statesLogLikelihood, THRESHOLD);
			
		}
		
		
		
	}

	@Test
	public void testCalculateStateLogLikelihoodFlat() throws Exception{
		
		StateLikelihood stateLikelihood = new ProbabilityStateLikelihood(); 
		double[] expecteds = new double[] { -4.4499971717799,
				-2.81959647584712, -2.22045226345552, -1.84839331479565,
				-1.57784235084435, -1.36512012588848, -1.18980694885304,
				-1.04069249808234, -0.910954992284636, -0.796132740880369,
				-0.693147180559945, -0.599784350155737, -0.514398666152773,
				-0.43573361212845, -0.362807998442184, -0.294840969650261,
				-0.231200924941943, -0.17136974738986, -0.114917146243339,
				-0.0614818641436684, -0.0117473304903028 };
		for (int i = 0; i < frequencies.length; i++) {
			double log = stateLikelihood.caluclateStateLogLikelihood(frequencies[i]);
			assertEquals("fail at "+i, expecteds[i], log, THRESHOLD);
		}
		
	}
	@Test
	public void testCalculateStateLogLikelihoodBetaMean() throws Exception{
		
		StateLikelihood stateLikelihood = new BetaMeanStateLikelihood();
		double[] expecteds = new double[] { -4.4627970966445,
				-4.45490108164499, -4.4088290544743, -4.35662071231919,
				-4.29972297213143, -4.23826264916909, -4.17195884206474,
				-4.10029323748722, -4.02253577268687, -3.93771569591255,
				-3.84455269254967, -3.74133935333329, -3.62574761789222,
				-3.49450146814089, -3.34279315588983, -3.16316066347848,
				-2.9430953103395, -2.65914011951639, -2.25862508309307,
				-1.57347309663271, 2.29615312974098 };
		for (int i = 0; i < 21; i++) {
			double log = stateLikelihood.caluclateStateLogLikelihood(frequencies[i]);
			assertEquals("fail at "+i, expecteds[i], log, THRESHOLD);
		}
	}
	
	@Test
	public void testCalculateStateLogLikelihoodBetaMode() throws Exception{
		
		StateLikelihood stateLikelihood = new BetaModeStateLikelihood();
		double[] expecteds = new double[] { -6.13013650123304,
				-2.26051027485935, -1.57535828839899, -1.17484325197567,
				-0.890888061152559, -0.670822708013584, -0.491190215602235,
				-0.339481903351173, -0.20823575359984, -0.0926440181587687,
				0.0105693210576128, 0.103732324420493, 0.18855240119481,
				0.266309865995157, 0.337975470572683, 0.404279277677028,
				0.465739600639366, 0.522637340827131, 0.574845682982234,
				0.620917710152927, 0.628813725152435 };
		for (int i = 0; i < 21; i++) {
			double log = stateLikelihood.caluclateStateLogLikelihood(frequencies[i]);
			assertEquals("fail at " + i, expecteds[i], log, THRESHOLD);
		}
	}
	

	@Test
	public void testCalculateStateLogLikelihoodGTest() throws Exception{
		
		StateLikelihood stateLikelihood = new GTestStateLikelihood();
		double[] expecteds = new double[] { -6.54540700719581,
				-6.08520556397393, -5.68744195480934, -5.31842484136728,
				-4.96806776933106, -4.63151576637113, -4.30581779805705,
				-3.98890554821967, -3.67915590009021, -3.37515321539537,
				-3.07552433951212, -2.77878684442851, -2.48316519890541,
				-2.18631510268386, -1.88484260502715, -1.57335518321388,
				-1.24232964427087, -0.872435902109154, -0.414878166309604,
				0.323417013652523, 1.18179430054454 };
		for (int i = 0; i < 21; i++) {
			double log = stateLikelihood.caluclateStateLogLikelihood(frequencies[i]);
			assertEquals("fail at " + i, expecteds[i], log, THRESHOLD);
		}
	}
	

	@Test
	public void testCalculateStateLogLikelihoodChisq() throws Exception{
		

		StateLikelihood stateLikelihood = new ChisqStateLikelihood();

		double[] expecteds = new double[] { -49.3168653997597,
				-44.8046134842287, -40.4312738560728, -36.2909399823289,
				-32.3832226986611, -28.7076563643503, -25.2636774106532,
				-22.0505947966424, -19.0675484079555, -16.3134489648306,
				-13.786888606292, -11.4860031036221, -9.40825045103395,
				-7.55003641111874, -5.90603939335547, -4.46788851746636,
				-3.22127082107978, -2.13850678309654, -1.15398942877333,
				-0.0294953974800052, 1.43811198038661 };
		for (int i = 0; i < 21; i++) {
			double log = stateLikelihood.caluclateStateLogLikelihood(frequencies[i]);
			assertEquals("fail at " + i, expecteds[i], log, THRESHOLD);
		}
	}
	
}
