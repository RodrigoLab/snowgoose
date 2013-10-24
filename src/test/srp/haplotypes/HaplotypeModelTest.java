package test.srp.haplotypes;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import org.apache.commons.lang3.StringUtils;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.Haplotype;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.HaplotypeModelUtils;
import srp.haplotypes.SwapInfo;
import srp.haplotypes.operator.SingleBaseOperator;
import srp.haplotypes.operator.MultiBasesOperator;
import srp.likelihood.LikelihoodUtils;
import test.TestUtils;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.util.Taxon;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;


public class HaplotypeModelTest {

	private static AlignmentMapping aMap;
	private static Alignment srpAlignment;
	
	private static Taxon[] expectedTaxons;
	private static String[] expectedSequences;
	private static ArrayList<Haplotype> expectedList;
//	private static HashMap<Character, Integer> charToState;

	
	private HaplotypeModel haplotypeModel;
	private HaplotypeModel haplotypeModelRandom;

	
	
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {

		String dir = System.getProperty("user.dir")+File.separatorChar+"unittest"+File.separator;
		srpAlignment = DataImporter.importAlignment(dir, "HaplotypeModelTest_10.fasta");
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
		expectedTaxons = new Taxon[10];
//		char[][] expectedMatrix = new char[matrixS.length][matrixS[0].length()];
		expectedList = new ArrayList<Haplotype>();
		for (int i = 0; i < expectedSequences.length; i++) {
			expectedTaxons[i] = new Taxon("r"+(i+1)+".1");
			Haplotype h = new Haplotype(expectedTaxons[i], expectedSequences[i]);
			expectedList.add(h);
		}
		
//		charToState = new HashMap<Character, Integer>();
//		charToState.put('A', 0);
//		charToState.put('C', 1);
//		charToState.put('G', 2);
//		charToState.put('T', 3);
//		charToState.put('-', 17);
//		charToState.put('*', 17);
//		charToState.put('.', 17);
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		haplotypeModel = new HaplotypeModel(aMap, srpAlignment);
		haplotypeModelRandom = new HaplotypeModel(aMap, 5);
	}

	@After
	public void tearDown() throws Exception {
	}

	
	
	
	@Test
	public void testInit() throws Exception {
		
		assertEquals(10, haplotypeModel.getHaplotypeCount());
		assertEquals(10, haplotypeModel.getSequenceCount());
		assertEquals(10, haplotypeModel.getTaxonCount());
		assertEquals(100, haplotypeModel.getHaplotypeLength());
		for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
			Taxon expected = new Taxon("r"+(i+1)+".1");
			assertEquals(expected, haplotypeModel.getTaxon(i));
		}
		

		assertEquals(5, haplotypeModelRandom.getHaplotypeCount());
		assertEquals(5, haplotypeModelRandom.getSequenceCount());
		assertEquals(5, haplotypeModelRandom.getTaxonCount());
		assertEquals(100, haplotypeModelRandom.getHaplotypeLength());
		for (int i = 0; i < haplotypeModelRandom.getHaplotypeCount(); i++) {
			Taxon expected = new Taxon(HaplotypeModel.TAXON_PREFIX+i);
			assertEquals(expected, haplotypeModelRandom.getTaxon(i));
		}
	}
	
	@Test
	public void testGetAlignment() throws Exception {
		
	
		String[] newHaplotypes = expectedSequences;

		for (int i = 0; i < newHaplotypes.length; i++) {
			assertEquals(haplotypeModel.getHaplotypeString(i), haplotypeModel.getAlignedSequenceString(i)) ;
			assertEquals(newHaplotypes[i], haplotypeModel.getHaplotypeString(i));
			for (int j = 0; j < newHaplotypes[i].length(); j++) {
				char expected = newHaplotypes[i].charAt(j); 
				assertEquals(expected,haplotypeModel.getHaplotypeCharAt(i,j));
			}
			
		}
	}
	@Test
	public void testReject() throws Exception {
		
		MCMCOperator op;
		op = new SingleBaseOperator(haplotypeModel, 0);
		
		for (int i = 0; i < 1000; i++) {
			op.operate();
			op.reject();
			haplotypeModel.reject();
		}
		for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
			assertEquals(haplotypeModel.getAlignedSequenceString(i), srpAlignment.getAlignedSequenceString(i));
		}
		
		op = new MultiBasesOperator(haplotypeModel, 10, null);
		for (int i = 0; i < 1000; i++) {
			op.operate();
			op.reject();
			haplotypeModel.reject();
		}
		for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
			assertEquals(haplotypeModel.getAlignedSequenceString(i), srpAlignment.getAlignedSequenceString(i));
		}
	
	}
	

	@Test
	public void testSetupAlignment() throws Exception {

		haplotypeModel = new HaplotypeModel(aMap, 5);
		for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
			assertNotEquals(haplotypeModel.getHaplotypeString(i), srpAlignment.getAlignedSequenceString(i));
			assertNotSame(haplotypeModel.getHaplotypeString(i), srpAlignment.getSequence(i));
			assertNotEquals(haplotypeModel.getHaplotypeString(i), expectedSequences[0] );
			assertNotSame(haplotypeModel.getHaplotypeString(i), expectedSequences[0]);
		
			assertNotSame(haplotypeModel.getHaplotypeString(i), srpAlignment.getSequence(i));
			assertNotEquals(haplotypeModel.getHaplotypeString(i), expectedSequences[0] );
			assertNotSame(haplotypeModel.getHaplotypeString(i), expectedSequences[0]);
		}	
	}

	@Test
	public void testHaplotypeGetNextBaseFrequency() throws Exception {


		Parameter frequency = new Parameter.Default(new double[]{0.1,0.25,0.60,0.05});

		int[] count = new int['T'+1]; 
		int ite = 1000;
    	for (int i = 0; i < ite; i++) {
    		int[] posChar = haplotypeModel.getNextBaseFrequency(frequency);
    		count[posChar[1]]++;
		}
    	TestUtils.assertExpectationRange(count['A'], frequency.getParameterValue(0)*ite, 0.05*ite);
    	TestUtils.assertExpectationRange(count['C'], frequency.getParameterValue(1)*ite, 0.05*ite);
    	TestUtils.assertExpectationRange(count['G'], frequency.getParameterValue(2)*ite, 0.05*ite);
    	TestUtils.assertExpectationRange(count['T'], frequency.getParameterValue(3)*ite, 0.05*ite);

    	
    	
    	double[] logFreq = new double[4];
    	for (int i = 0; i < logFreq.length; i++) {
			logFreq[i] = Math.log(frequency.getParameterValue(i));
		}
    	for (int i = 0; i < logFreq.length; i++) {
			for (int j = 0; j < logFreq.length; j++) {
				double expected = logFreq[i]-logFreq[j];
				assertEquals(expected, haplotypeModel.getLogqFrequencyStates(i,j), 0);
			}
		}
    	
    	frequency = new Parameter.Default(new double[]{0.1, 0.2, 0.3, 0.4});
    	haplotypeModel.getNextBaseFrequency(frequency);
     	
    	for (int i = 0; i < logFreq.length; i++) {
			logFreq[i] = Math.log(frequency.getParameterValue(i));
		}
    	for (int i = 0; i < logFreq.length; i++) {
			for (int j = 0; j < logFreq.length; j++) {
				double expected = logFreq[i]-logFreq[j];
				assertEquals(expected, haplotypeModel.getLogqFrequencyStates(i,j), 0);
			}
		}
    	count = new int['T'+1]; 
    	for (int i = 0; i < ite; i++) {
    		int[] posChar = haplotypeModel.getNextBaseFrequency(frequency);
    		count[posChar[1]]++;
		}
    	TestUtils.assertExpectationRange(count['A'], frequency.getParameterValue(0)*ite, 0.05*ite);
    	TestUtils.assertExpectationRange(count['C'], frequency.getParameterValue(1)*ite, 0.05*ite);
    	TestUtils.assertExpectationRange(count['G'], frequency.getParameterValue(2)*ite, 0.05*ite);
    	TestUtils.assertExpectationRange(count['T'], frequency.getParameterValue(3)*ite, 0.05*ite);

    	
    	
	}
	
	
	@Test
	public void testCalculateSpsSingle() {
		String[] seqs = new String[]{
				"AAAAACCCCCGGGGGTTTT",
				"AAAACCCCCGGGGGTTTTA"};
		Alignment alignment = AlignmentUtils.createAlignment(seqs);
		HaplotypeModel haplotype = new HaplotypeModel(AlignmentUtils.createAlignmentMapping(seqs), alignment);
		int expected = 4;
		assertEquals(expected, haplotype.calculateSPS());

		seqs = new String[]{
				"AAAAACCCCCGGGGGTTTT",
				"AAAACCCCCGGGGGTTTTA",
				"AAACCCCCGGGGGTTTTAA"
		};
		alignment = AlignmentUtils.createAlignment(seqs);
		haplotype = new HaplotypeModel(AlignmentUtils.createAlignmentMapping(seqs), alignment);
		expected = 4+4+8;
		assertEquals(expected, haplotype.calculateSPS());
		

		seqs = new String[]{
				"ACGTTATTTT",
				"ACGTTATTTT",
				"ACGTTATTTT",
		};
		alignment = AlignmentUtils.createAlignment(seqs);
		haplotype = new HaplotypeModel(AlignmentUtils.createAlignmentMapping(seqs), alignment);
		expected = 0;
		assertEquals(expected, haplotype.calculateSPS());
		
		seqs = new String[]{
				"ACGTTATTTT",
				"CGTATTATTT",
				"GTACTTTATT"
		};
		alignment = AlignmentUtils.createAlignment(seqs);
		haplotype = new HaplotypeModel(AlignmentUtils.createAlignmentMapping(seqs), alignment);
		expected = 3+3+3+3+2+2+2;
		assertEquals(expected, haplotype.calculateSPS());
	
	}

	@Test
	public void testCalculateSpsDoubleSame() {
	
		String[] seq1 = new String[]{
				"AACCTTGG",
				"ACCTTGGA",
				"CCTTGGAA"
		};
		String[] seq2 = new String[]{
				"AACCTTGG",
				"ACCTTGGA",
				"CCTTGGAA"
		};
		HaplotypeModel h1 = AlignmentUtils.createHaplotypeModel(seq1);
		HaplotypeModel h2 = AlignmentUtils.createHaplotypeModel(seq2);
		int sps = HaplotypeModelUtils.calculeteSPS(h1, h2);
		int expected = 4+8 + 4+4 + 4+8;
		assertEquals(expected, sps, 0);
		
		int[][] spsArray = HaplotypeModelUtils.calculeteSPSArray(h1, h2);
		int[][] expectedArray = new int[][]{
					{ 0, 4, 8 }, 
					{ 4, 0, 4 },
					{ 8, 4, 0 }		
				};
		
		seq2 = new String[]{
					"AACTCTGA", //3+3+6
					"AAGTCCTA", //6+5+6
					"AACAATGA"	//3+4+7
				};
		h2 = AlignmentUtils.createHaplotypeModel(seq2);
		sps = HaplotypeModelUtils.calculeteSPS(h1, h2);
		expected = 43;
		assertEquals(expected, sps);

		spsArray = HaplotypeModelUtils.calculeteSPSArray(h1, h2);
		expectedArray = new int[][] { 
						{ 3, 6, 3 }, 
						{ 3, 5, 4 }, 
						{ 6, 6, 7 } 
				};
		assertArrayEquals(expectedArray, spsArray);
	}

	
	
}
