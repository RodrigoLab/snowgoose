package test.srp.haplotypes;

import static org.junit.Assert.assertEquals;


import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.ShortRead;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

public class AlignmentMappingTest {

	private static Alignment shortReads;

	private AlignmentMapping aMap;

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		
		String dir = System.getProperty("user.dir")+File.separatorChar+"unittest"+File.separator;
		String srpFileName = "align_test_4.fasta";
		shortReads = DataImporter.importAlignment(dir, srpFileName);
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
