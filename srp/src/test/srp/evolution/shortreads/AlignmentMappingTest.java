package test.srp.evolution.shortreads;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.ArrayList;

import org.apache.commons.math3.stat.StatUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.evolution.haplotypes.AlignmentUtils;
import srp.evolution.shortreads.AlignmentMapping;
import srp.evolution.shortreads.ShortRead;
import test.TestUtils;
import dr.evolution.alignment.Alignment;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

public class AlignmentMappingTest {

	private static Alignment shortReads;

	private AlignmentMapping aMap;

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		
		String dir = TestUtils.getUnittestDir();
		String srpFileName = "AlignmentMappingTest_4_srp.fasta";
		shortReads = DataImporter.importShortReads(dir, srpFileName);
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {

		aMap = new AlignmentMapping(shortReads);
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testAlignmentMapping() {

		
		assertEquals(24, aMap.getLength() );
		assertEquals(4, aMap.getSrpCount());
		
		String[] seqs = new String[] {
			"AA******AA********A.....",
			"..CC****CCCC******CC....",
			"....GG**GGGGGG****GGG...",
			"......TTTTTTTTTT**TTTT.."};
		String[] taxons = new String[]{
			"srpA","srpC","srpG","srpT"
		};
		
		for (int i = 0; i < aMap.getSrpCount(); i++) {
			Sequence s = new Sequence(new Taxon(taxons[i]), seqs[i]);
			ShortRead srp = new ShortRead(s);
			assertEquals(srp.getStart(), aMap.getSrpStart(i));
			assertEquals(srp.getEnd(), aMap.getSrpEnd(i));
			assertEquals(srp.getFragmentSrp(), aMap.getSrpFragment(i));
			assertEquals(srp.getFullSrp(), aMap.getSrpFull(i));
			assertEquals(srp.getName(), aMap.getSrpName(i));
			assertEquals(srp.getLength(), aMap.getSrpLength(i));
		}
		
		
	}
	
	@Test
		public void testGetNextBase() {
	
	
			String[] seqs = new String[] {
				"A***AAA...",
				".C**AAA...",
				"..G*AAA...",
				"...TTGC..."};
			
			aMap = new AlignmentMapping(AlignmentUtils.createAlignment(seqs));
					
			for (int i = 0; i < 1e5; i++) {
				int[] posChar = aMap.getNextBase();
				switch (posChar[0]) {
				case 0:
					assertEquals(posChar[1], 'A');
					break;
				case 1:
					assertTrue( posChar[1]=='C' || posChar[1]=='-');
					break;
				case 2:
					assertTrue( posChar[1]=='G' || posChar[1]=='-');
					break;
				case 3:
					assertTrue( posChar[1]=='T' || posChar[1]=='-');
					break;
				case 4:
					assertTrue( posChar[1]=='A' || posChar[1]=='T');
					break;
				case 7:
					assertEquals(posChar[1], '-');
					break;
				default:
					break;
				}
			}
			
			seqs = new String[] {
					"AAA",
					"AAA",
					"AAA",
					"AAA",
					"AAA",
					"AAA",
					"AAA",
					"AAA",
					"AAA",
					"AAA",
					"CGT"};
	
			aMap = new AlignmentMapping(AlignmentUtils.createAlignment(seqs));
			int[][] count = new int[3][85];
			for (int i = 0; i < 3e5; i++) {
				int[] posChar = aMap.getNextBase();
				switch (posChar[0]) {
					case 0:
						count[0][posChar[1]]++;
						break;
					case 1:
						count[1][posChar[1]]++;
						break;
					case 2:
						count[2][posChar[1]]++;
						break;
					default:
						break;
				}
			}
			double expectedRatio = 10;
			double ratio = (double)count[0]['A'] / (double)count[0]['C'];
			TestUtils.assertExpectationRange(ratio, expectedRatio , 0.25);
			ratio = (double)count[1]['A'] / (double)count[1]['G'];
			TestUtils.assertExpectationRange(ratio, expectedRatio , 0.25);
			ratio = (double)count[2]['A'] / (double)count[2]['T'];
			TestUtils.assertExpectationRange(ratio, expectedRatio , 0.25);
		}



	@Test
	public void testGetNextBaseUniform() {

	
		String[] seqs = new String[] {
				"AAAAA.",
				"AAAAA.",
				"AAAAA.",
				"AAAAA.",
				"AAAAA.",
				"AAAAA.",
				"AAAAA.",
				"AAAAA.",
				"AAACA.",
				"AAAGA.",
				"CGTT*."};

		aMap = new AlignmentMapping(AlignmentUtils.createAlignment(seqs));
		double[][] count = new double[6][85];
		for (int i = 0; i < 5e5; i++) {
			int[] posChar = aMap.getNextBaseUniform();
			count[posChar[0]][posChar[1]]++;
			
		}
		double expectedRatio = 1;
		double ratio = count[0]['A'] / count[0]['C'];
		TestUtils.assertExpectationRange(ratio, expectedRatio , 0.1);
		ratio = count[1]['A'] / count[1]['G'];
		TestUtils.assertExpectationRange(ratio, expectedRatio , 0.1);
		ratio = count[2]['A'] / count[2]['T'];
		TestUtils.assertExpectationRange(ratio, expectedRatio , 0.1);
		
		ratio = count[3]['A'] / count[3]['C'];
		TestUtils.assertExpectationRange(ratio, expectedRatio , 0.1);
		ratio = count[3]['A'] / count[3]['G'];
		TestUtils.assertExpectationRange(ratio, expectedRatio , 0.1);
		ratio = count[3]['A'] / count[3]['T'];
		TestUtils.assertExpectationRange(ratio, expectedRatio , 0.1);
		
		ratio = count[4]['A'] / count[4]['-'];
		assertEquals(StatUtils.sum(count[4]), count[4]['A'], 0);
		assertEquals(count[4]['-'], 0, 0);
		
		ratio = count[4]['A'] / count[5]['-'];
		assertEquals(StatUtils.sum(count[5]), count[5]['-'], 0);

		
	}
	

	@Test
	public void testGetNextBaseEmpirical() {

		String[] seqs = new String[] {
				"AAATT.",
				"AAAAA.",
				"AAAAA.",
				"AAAAA.",
				"AAACC.",
				"AAACC.",
				"ACCAA.",
				"ACCAA.",
				"CAAAA.",
				"CAAAA.",
				"TTTAA."};
//A:C:T = 8:2:1
		aMap = new AlignmentMapping(AlignmentUtils.createAlignment(seqs));
		double[][] count = new double[6][85];
		for (int i = 0; i < 1e6; i++) {
			int[] posChar = aMap.getNextBaseEmpirical();//FIXME
			count[posChar[0]][posChar[1]]++;
		}
		
		double ratio = 1;
		for (int i = 0; i < count.length; i++) {
			ratio = count[i]['C'] / count[0]['A'];
			TestUtils.assertExpectationRange(ratio, 0.25, 0.01);
			ratio = count[i]['G'] / count[1]['A'];
			TestUtils.assertExpectationRange(ratio, 0, 1e-8);
			ratio = count[i]['T'] / count[2]['A'];
			TestUtils.assertExpectationRange(ratio, 0.125, 0.01);
		}
		
	}


	@Test
	public void testMapNameToID() {
		String[] taxons = new String[]{
				"srpA","srpC","srpG","srpT"
		};
		for (int i = 0; i < taxons.length; i++) {
			Integer actual= aMap.mapNameToID(taxons[i]);
			assertEquals(Integer.valueOf(i), actual);
		}
	}
	
	@Test
	public void testMapToSrp(){
//		String[] seqs = new String[] {
//				"AA******AA********A.....",
//				"..CC****CCCC******CC....",
//				"....GG**GGGGGG****GGG...",
//				"......TTTTTTTTTT**TTTT.."};
		ArrayList<Integer> expected = new ArrayList<Integer>();
		for (int i = 0; i < 4; i++) {
			expected.add(i);
			assertEquals(expected, aMap.getMapToSrp(2*i));
			assertEquals(expected, aMap.getMapToSrp(2*i+1));
		}
		for (int i = 0; i < 4; i++) {
			assertEquals(expected, aMap.getMapToSrp(18+i));
			expected.remove(0);
		}
		assertEquals(expected, aMap.getMapToSrp(22));
	}
}
