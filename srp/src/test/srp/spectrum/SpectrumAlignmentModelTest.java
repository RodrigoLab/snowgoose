package test.srp.spectrum;

import static org.junit.Assert.*;

import java.io.File;
import java.util.ArrayList;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import dr.evolution.alignment.Alignment;
import dr.evolution.util.Taxon;

import srp.core.DataImporter;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.Haplotype;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;

public class SpectrumAlignmentModelTest {
	private AlignmentMapping aMap;
	private String[] expectedSequences;

	@Before
	public void createSpectrumAlignmentModel() throws Exception {

		String dir = System.getProperty("user.dir")+File.separatorChar+"unittest"+File.separator;
		Alignment srpAlignment = DataImporter.importShortReads(dir, "HaplotypeModelTest_10_srp.fasta");
		aMap = new AlignmentMapping(srpAlignment);
		

		expectedSequences = new String[]{ 
				".............................................AAGCTCCAACACA*GAGAGCAACTCGTAGGGCGTTGTTCAA..............",
				"..............CCGCGCC*TGTCAGGGCTCAAAACCC*GA*AAAGC...................................................",
				"..............................CTCAAAATCCTGAGAAAGCTCCAACGCACGAGAGCAACT...............................",
				".................CGCC*TGTCAGGGCTCAAAATCCTGAGAAAGCTCC................................................",
				".........................TAGGGCTCAAAATCCTGAGAAAGCTCCAACACACGAG......................................",
				"CCAGG...............................................................................................",
				"...............................................GCTCCAACACACGAGAGC*ACTCGTAGGGCGTTGTTCAATTC...........",
				".......GTCCAGTCCGCGCCTTGTCAGGGCTCG..................................................................",
				"........................................................CACGAGAGCAACACGTA*GGCGTTG*TCAATTCTAGTTCT....",
				"..........CAGTCCGCGCCTTGTCA*GGCTCAAAAT*CTG..........................................................",
		};
//		expectedTaxons = new Taxon[10];
////		char[][] expectedMatrix = new char[matrixS.length][matrixS[0].length()];
//		expectedList = new ArrayList<Haplotype>();
//		for (int i = 0; i < expectedSequences.length; i++) {
//			expectedTaxons[i] = new Taxon("r"+(i+1)+".1");
//			Haplotype h = new Haplotype(expectedTaxons[i], expectedSequences[i]);
//			expectedList.add(h);
//		}
		
	}

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
	public void testConstructor() throws Exception {
		double[] expectedFreq = new double[]{0.25, 0.25, 0.25, 0.25};
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 5);
		assertEquals(100, spectrumModel.getSpectrumLength());
		assertEquals(5, spectrumModel.getSpectrumCount());
		for (int i = 0; i < spectrumModel.getSpectrumCount(); i++) {
			Spectrum spectrum = spectrumModel.getSpectrum(i);
			assertArrayEquals(expectedFreq , spectrum.getFrequenciesAt(0), 0);
		}
		
	}
	

	@Test
	public void testStoreRestore() throws Exception {
		
	
	}
}
