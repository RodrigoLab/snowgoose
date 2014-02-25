package test.srp.spectrum;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.haplotypes.AlignmentMapping;
import srp.spectrum.SpectraParameter;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumAlignmentUtils;
import dr.evolution.alignment.Alignment;

public class SpectrumAlignmentUtilsTest {

	private static final double MIN_FREQ = SpectraParameter.MIN_FREQ;
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}


	private double THRESHOLD = 1e-6;

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	
	@Test
	public void testDist() throws Exception {
		String dataDir;

		dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/";

		String trueHaplotypeFile = "H7_fullHaplotype.fasta";
		String logSpectrumPath = dataDir+"H7_haplotypepartial_01.fasta";

		DataImporter dataImporter = new DataImporter(dataDir);
		Alignment trueAlignment = dataImporter.importAlignment(trueHaplotypeFile);
		AlignmentMapping aMap = new AlignmentMapping(trueAlignment);

		SpectrumAlignmentModel spectrumModel = dataImporter.importPartialSpectrumFile(aMap, logSpectrumPath );
//			
		SpectrumAlignmentUtils.compareSpectrumToTrueAlignment(spectrumModel, trueAlignment);

	}
	
	
	@Test
	public void testDist01() throws Exception {

		String dataDir;

		dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/";
		String trueHaplotypeFile = "spectrumDist_true.fasta";
//		>hap_0	GTGCA
//		>hap_1	ATGTA
//		>hap_2	GTGCA
//		>hap_3	GTGCA
		String partialSpectrumFile = "spectrumDist_partial.fasta"; //GTGCA
									  
		
		DataImporter dataImporter = new DataImporter(dataDir);
		Alignment trueAlignment = dataImporter.importAlignment(trueHaplotypeFile);
		
		AlignmentMapping aMap = new AlignmentMapping(trueAlignment);
		SpectrumAlignmentModel spectrumModel = dataImporter.importPartialSpectrumFile(aMap, partialSpectrumFile );

		double[][] allDelta = SpectrumAlignmentUtils.compareSpectrumToTrueAlignment(spectrumModel, trueAlignment);
//		Spectrum1
//		>hap_0	GTGCA MMMMM
//		>hap_1	ATGTA XMMXM
//		>hap_2	AAACA XXXMM
//		>hap_3	TATGC XXXXX
		
		double match = Math.sqrt(0.3*0.3 + 3*0.1*0.1);
		double mismatch = Math.sqrt(0.9*0.9 + 0.7*0.7 + 2*0.1*0.1);
		double[][] expectedDist = new double[][]{
				{match*5, match*3+mismatch*2, match*2+mismatch*3, mismatch*5 },
				{mismatch, match, mismatch, match, },  
				{match, match, match, match, },  
				{mismatch, mismatch, mismatch,  mismatch*5 }  
				
		};
		double[][] allDeltas = allDelta;
		assertArrayEquals(expectedDist[0], allDelta[0], THRESHOLD );
		System.out.println(SpectrumAlignmentUtils.formatter(allDeltas));
		System.out.println(SpectrumAlignmentUtils.formatter(expectedDist));
		
		
	}

}
