package test.srp.evolution.spectrum;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.evolution.spectrum.SpectraParameter;
import srp.evolution.spectrum.SpectraParameter.SpectraType;
import srp.evolution.spectrum.Spectrum;
import srp.evolution.spectrum.SpectrumAlignmentModel;
//import srp.spectrum.AbstractSpectra;
//import srp.spectrum.AbstractSpectrum;
import dr.evolution.alignment.Alignment;

public class SpectrumAlignmentModelTest {
	private int spectrumLength;
	private String[] expectedSequences;

	@Before
	public void createSpectrumAlignmentModel() throws Exception {

		String dir = System.getProperty("user.dir")+File.separatorChar+"unittest"+File.separator;
		Alignment srpAlignment = DataImporter.importShortReads(dir, "HaplotypeModelTest_10_srp.fasta");

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
		spectrumLength = expectedSequences[0].length();
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
	public void testConstructorEqual() throws Exception {
		double[] expectedFreq = new double[]{0.25, 0.25, 0.25, 0.25};
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, 5, SpectraType.EQUAL);
		assertEquals(100, spectrumModel.getSpectrumLength());
		assertEquals(5, spectrumModel.getSpectrumCount());
		for (int i = 0; i < spectrumModel.getSpectrumCount(); i++) {
			Spectrum spectrum = spectrumModel.getSpectrum(i);
			assertArrayEquals(expectedFreq , spectrum.getFrequenciesAt(0), 0);
			assertArrayEquals(expectedFreq , spectrum.getFrequenciesAt(1), 0);
			assertArrayEquals(expectedFreq , spectrum.getFrequenciesAt(2), 0);
			assertArrayEquals(expectedFreq , spectrum.getFrequenciesAt(3), 0);
		}
		
	}
	
	@Test
	public void testConstructorZeroOne() throws Exception {
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, 5, SpectraType.DOMINANT);
		assertEquals(100, spectrumModel.getSpectrumLength());
		assertEquals(5, spectrumModel.getSpectrumCount());
		for (int i = 0; i < spectrumModel.getSpectrumCount(); i++) {
			Spectrum spectrum = spectrumModel.getSpectrum(i);
			for (int j = 0; j < spectrum.getLength(); j++) {
				SpectraParameter spectra = spectrum.getSpectra(j);
				int count0 = 0;
				int count1 = 0;
				for (int k = 0; k < spectra.getDimension(); k++) {
					double freq = spectra.getFrequency(k);
					if(freq==SpectraParameter.INIT_SMALL){
						count0++;
					}
					else if(freq==SpectraParameter.INIT_LARGE){
						count1++;
					}
				}
				assertEquals(3, count0);
				assertEquals(1, count1);
			}
		}
		
	}


	@Test
	public void testConstructorRandom() throws Exception {
		
		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(spectrumLength, 5, SpectraType.RANDOM);
		assertEquals(100, spectrumModel.getSpectrumLength());
		assertEquals(5, spectrumModel.getSpectrumCount());
		for (int i = 0; i < spectrumModel.getSpectrumCount(); i++) {
			Spectrum spectrum = spectrumModel.getSpectrum(i);
			for (int j = 0; j < spectrum.getLength(); j++) {
				SpectraParameter spectra = spectrum.getSpectra(j);
				double[] frequencies = spectra.getFrequencies();
				int count = 0;
				for (int k = 0; k < frequencies.length; k++) {
					for (int k2 = k+1; k2 < frequencies.length; k2++) {
						if(frequencies[k]!=frequencies[k2]){
							count++;
						}
					}
					
				}
				assertEquals(Arrays.toString(frequencies), 6, count);
			}
		}
		
	}
	
	
//	@Test
//	public void testGetSpectrum() throws Exception {
//	
//		
//	}
	
	@Test
	public void testStoreRestore() throws Exception {
		//TODO implement testStoreRestore()
	}

}
