package test.srp.evolution.haplotypes;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotSame;

import java.io.File;
import java.util.ArrayList;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.evolution.haplotypes.AlignmentUtils;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import srp.evolution.haplotypes.SPSDist;
import srp.evolution.shortreads.ShortReadMapping;
//import srp.evolution.haplotypes.old.OldHaplotype;
//import srp.evolution.haplotypes.old.OldHaplotypeModel;
//import srp.evolution.haplotypes.old.OldHaplotypeModelUtils;
//import srp.evolution.shortreads.AlignmentMapping;
import srp.operator.haplotypes.BaseSingleOperator;
import srp.operator.haplotypes.BasesMultiOperator;
import test.TestUtils;
import dr.evolution.alignment.Alignment;
import dr.evolution.util.Taxon;
import dr.inference.operators.MCMCOperator;
//import srp.operator.haplotypes.old.BaseSingleOperator;
//import srp.operator.haplotypes.old.BasesMultiOperator;


public class HaplotypeModelTest {

	private static ShortReadMapping aMap;
	private static Alignment srpAlignment;
	
	private static Taxon[] expectedTaxons;
	private static String[] expectedSequences;
	private static ArrayList<Haplotype> expectedList;
//	private static HashMap<Character, Integer> charToState;

	
	private HaplotypeModel haplotypeModel;
	private HaplotypeModel haplotypeModelRandom;

	
	
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {

		String dir = TestUtils.getUnittestDir();
		srpAlignment = DataImporter.importShortReads(dir, "HaplotypeModelTest_10_srp.fasta");
		aMap = new ShortReadMapping(srpAlignment);
		

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
		expectedTaxons = new Taxon[expectedSequences.length];
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
		haplotypeModel = new HaplotypeModel(srpAlignment);
		haplotypeModelRandom = new HaplotypeModel(5, 100);
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
		Haplotype haplotype = haplotypeModel.getHaplotype(0);
		System.out.println(haplotype.getChar(0));
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
		op = new BaseSingleOperator(haplotypeModel);
		
		for (int i = 0; i < 1000; i++) {
			haplotypeModel.storeModelState();
			op.operate();
			op.reject();
			haplotypeModel.restoreModelState();
		}
		for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
			assertEquals(haplotypeModel.getAlignedSequenceString(i), srpAlignment.getAlignedSequenceString(i));
		}
		
		op = new BasesMultiOperator(haplotypeModel, 10, null);
		for (int i = 0; i < 1000; i++) {
			haplotypeModel.storeModelState();
			op.operate();
			op.reject();
			haplotypeModel.restoreModelState();
		}
		for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
			assertEquals(haplotypeModel.getAlignedSequenceString(i), srpAlignment.getAlignedSequenceString(i));
		}
	
	}
	

	@Test
	public void testSetupAlignment() throws Exception {

		haplotypeModel = new HaplotypeModel(5, aMap.getLength());
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
	public void testCalculateSpsSingle() {
		String[] seqs = new String[]{
				"AAAAACCCCCGGGGGTTTT",
				"AAAACCCCCGGGGGTTTTA"};
		Alignment alignment = AlignmentUtils.createAlignment(seqs);
		HaplotypeModel haplotype = new HaplotypeModel(alignment);
		int expected = 4*2;
		assertEquals(expected, haplotype.calculateSPS());

		seqs = new String[]{
				"AAAAACCCCCGGGGGTTTT",
				"AAAACCCCCGGGGGTTTTA",
				"AAACCCCCGGGGGTTTTAA"
		};
		alignment = AlignmentUtils.createAlignment(seqs);
		haplotype = new HaplotypeModel( alignment);
		expected = 2*(4+4+8);
		assertEquals(expected, haplotype.calculateSPS());
		

		seqs = new String[]{
				"ACGTTATTTT",
				"ACGTTATTTT",
				"ACGTTATTTT",
		};
		alignment = AlignmentUtils.createAlignment(seqs);
		haplotype = new HaplotypeModel(alignment);
		expected = 0;
		assertEquals(expected, haplotype.calculateSPS());
		
		seqs = new String[]{
				"ACGTTATTTT",
				"CGTATTATTT",
				"GTACTTTATT"
		};
		alignment = AlignmentUtils.createAlignment(seqs);
		haplotype = new HaplotypeModel(alignment);
		expected = 2*(3+3+3+3+2+2+2);
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
		HaplotypeModel h1 = new HaplotypeModel(seq1);
		HaplotypeModel h2 = new HaplotypeModel(seq2);
		int sps = SPSDist.calculeteSPS(h1, h2);
		int expected = 4+8 + 4+4 + 4+8;
		assertEquals(expected, sps, 0);
		
		int[][] spsArray = SPSDist.calculeteSPSArray(h1, h2);
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
		h2 = new HaplotypeModel(seq2);
		sps = SPSDist.calculeteSPS(h1, h2);
		expected = 43;
		assertEquals(expected, sps);

		spsArray = SPSDist.calculeteSPSArray(h1, h2);
		expectedArray = new int[][] { 
						{ 3, 6, 3 }, 
						{ 3, 5, 4 }, 
						{ 6, 6, 7 } 
				};
		assertArrayEquals(expectedArray, spsArray);
	}

	
	
}
