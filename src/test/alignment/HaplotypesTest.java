package test.alignment;

import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import alignment.AlignmentMapping;
import alignment.AlignmentUtils;
import alignment.Haplotypes;
import core.DataImporter;
import dr.evolution.alignment.Alignment;

public class HaplotypesTest {

	private static AlignmentMapping aMap;

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		
		String dir = "/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/";
		String srpFileName = "1110_10_align_test.fasta";
		DataImporter dataImporter = new DataImporter(dir);

		Alignment shortReads = dataImporter.importAlignment(srpFileName);
//		aMap = new AlignmentMapping(shortReads);
//		refSeq = dataImporter.importRefSeq(refFileName);
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
	public void testCalculateSpsSingle() {
		String[] seqs = new String[]{
				"AAAAACCCCCGGGGGTTTT",
				"AAAACCCCCGGGGGTTTTA"};
		Alignment alignment = AlignmentUtils.createAlignment(seqs);
		Haplotypes haplotype = new Haplotypes(AlignmentUtils.createAlignmentMapping(seqs), alignment);
		double expected = 4;
		assertEquals(expected, haplotype.calculateSPS(), 0);

		seqs = new String[]{
				"AAAAACCCCCGGGGGTTTT",
				"AAAACCCCCGGGGGTTTTA",
				"AAACCCCCGGGGGTTTTAA"
		};
		alignment = AlignmentUtils.createAlignment(seqs);
		haplotype = new Haplotypes(AlignmentUtils.createAlignmentMapping(seqs), alignment);
		expected = 4+4+8;
		assertEquals(expected, haplotype.calculateSPS(), 0);
		

		seqs = new String[]{
				"ACGTTATTTT",
				"ACGTTATTTT",
				"ACGTTATTTT",
		};
		alignment = AlignmentUtils.createAlignment(seqs);
		haplotype = new Haplotypes(AlignmentUtils.createAlignmentMapping(seqs), alignment);
		expected = 0;
		assertEquals(expected, haplotype.calculateSPS(), 0);
		
		seqs = new String[]{
				"ACGTTATTTT",
				"CGTATTATTT",
				"GTACTTTATT"
		};
		alignment = AlignmentUtils.createAlignment(seqs);
		haplotype = new Haplotypes(AlignmentUtils.createAlignmentMapping(seqs), alignment);
		expected = 3+3+3+3+2+2+2;
		assertEquals(expected, haplotype.calculateSPS(), 0);
	
	}

	@Test
	public void testCalculateSpsDouble() {
	
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
		Haplotypes h1 = AlignmentUtils.createHaplotypes(seq1);
		Haplotypes h2 = AlignmentUtils.createHaplotypes(seq2);
		double sps = Haplotypes.calculeteSPS(h1, h2);
		double expected = 4+8 + 4+4 + 4+8;
		assertEquals(expected, sps, 0);
		
		seq2 = new String[]{
				"AACTCTGA", //3+3+6
				"AAGTCCTA", //6+5+6
				"AACAATGA"	//3+4+7
		};
		h2 = AlignmentUtils.createHaplotypes(seq2);
		sps = Haplotypes.calculeteSPS(h1, h2);
		expected = 43;
		assertEquals(expected, sps, 0);
		
	}

}
