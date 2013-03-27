package test.haplotypes;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.hamcrest.CoreMatchers;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import srp.core.DataImporter;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.Haplotype;
import srp.haplotypes.HaplotypeModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxon;


public class HaplotypeModelTest {

	private static AlignmentMapping aMap;
	private static Alignment srpAlignment;
	
	private static Taxon[] expectedTaxons;
	private static String[] expectedSequences;
	private static ArrayList<Haplotype> expectedList;
	private static HashMap<Character, Integer> charToState;
	
	private HaplotypeModel haplotypeModel;
	private HaplotypeModel haplotypeModelRandom;

	
	
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {

		String dir = System.getProperty("user.dir")+File.separatorChar+"unittest"+File.separator;
		srpAlignment = DataImporter.importAlignment(dir, "HaplotypesTest_10.fasta");
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
		
		charToState = new HashMap<Character, Integer>();
		charToState.put('A', 0);
		charToState.put('C', 1);
		charToState.put('G', 2);
		charToState.put('T', 3);
		charToState.put('-', 17);
		charToState.put('*', 17);
		charToState.put('.', 17);
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
	public void testHaplotypeExt() throws Exception {
		
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
		
	
		String[] newHaplotypes = new String[haplotypeModel.getHaplotypeCount()];
		for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
			haplotypeModel.randomSeq(i);
			newHaplotypes[i] = String.valueOf( haplotypeModel.getCharMatrix()[i] );
			assertEquals(newHaplotypes[i], haplotypeModel.getHaplotypeString(i));
		}

		Alignment expectedAlignment = haplotypeModel.getAlignment();
		for (int i = 0; i < newHaplotypes.length; i++) {
			assertEquals(newHaplotypes[i], expectedAlignment.getAlignedSequenceString(i)) ;
		}
	}
	@Test
	public void testReject() throws Exception {
		int seqIndex = 0;
		String seq = haplotypeModel.getHaplotypeString(seqIndex);
		for (int i = 0; i < 1000; i++) {
			haplotypeModel.swapBase(seqIndex);
			haplotypeModel.reject();
		}	
		String newSeq = haplotypeModel.getHaplotypeString(seqIndex);
		assertEquals(seq, newSeq);
		
		for (int i = 0; i < 1000; i++) {
			haplotypeModel.swapBase();
			haplotypeModel.reject();
		}
		for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
			assertEquals(haplotypeModel.getAlignedSequenceString(i), srpAlignment.getAlignedSequenceString(i));
		}
	
	}
	@Test
	public void testSwapBase() throws Exception {

		int seqIndex = 0;
		String seq = haplotypeModel.getHaplotypeString(seqIndex);
		for (int i = 0; i < 1000; i++) {
			haplotypeModel.swapBase(seqIndex);
		}	
		String newSeq = haplotypeModel.getHaplotypeString(seqIndex);
		
		assertNotEquals(seq, newSeq);
		assertEquals(6, StringUtils.countMatches(newSeq, "-")) ;

		seqIndex =1;
		for (int i = 0; i < 10; i++) {
			haplotypeModel.swapBase(seqIndex, i);
			assertEquals(newSeq.charAt(i),haplotypeModel.getHaplotypeString(seqIndex).charAt(i));
		}
		for (int i = 0; i < 1000; i++) {
			haplotypeModel.swapBase();
		}
		for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
			assertNotEquals(haplotypeModel.getAlignedSequenceString(i), srpAlignment.getAlignedSequenceString(i));
		}
	}
	@Test
	public void testSetupAlignment() throws Exception {

		String seq = haplotypeModel.getHaplotypeString(0);
		for (int i = 0; i < 1000; i++) {
			haplotypeModel.swapBase(0);
		}	
		String newSeq = haplotypeModel.getHaplotypeString(0);
		Alignment alignment = haplotypeModel.getAlignment();
		
		assertNotEquals(newSeq, srpAlignment.getAlignedSequenceString(0));
		assertNotSame(alignment.getSequence(0), srpAlignment.getSequence(0));
		assertNotSame(alignment, srpAlignment);
		
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
		int sps = HaplotypeModel.calculeteSPS(h1, h2);
		int expected = 4+8 + 4+4 + 4+8;
		assertEquals(expected, sps, 0);
		
		int[][] spsArray = HaplotypeModel.calculeteSPSArray(h1, h2);
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
		sps = HaplotypeModel.calculeteSPS(h1, h2);
		expected = 43;
		assertEquals(expected, sps);

		spsArray = HaplotypeModel.calculeteSPSArray(h1, h2);
		expectedArray = new int[][] { 
						{ 3, 6, 3 }, 
						{ 3, 5, 4 }, 
						{ 6, 6, 7 } 
				};
		assertArrayEquals(expectedArray, spsArray);
	}

	
	
}
