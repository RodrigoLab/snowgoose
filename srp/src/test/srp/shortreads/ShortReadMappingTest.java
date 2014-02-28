package test.srp.shortreads;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.ArrayList;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.shortreads.ShortRead;
import srp.shortreads.ShortReadMapping;
import dr.evolution.alignment.Alignment;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

public class ShortReadMappingTest {

	private static Alignment shortReads;

	private ShortReadMapping srpMap;

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		
		String dir = System.getProperty("user.dir")+File.separatorChar+"unittest"+File.separator;
		String srpFileName = "ShortReadMappingTest_4_srp.fasta";
		shortReads = DataImporter.importShortReads(dir, srpFileName);
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {

		srpMap = new ShortReadMapping(shortReads);
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testShortReadMapping() {

		
		assertEquals(24, srpMap.getLength() );
		assertEquals(4, srpMap.getSrpCount());
		
		String[] seqs = new String[] {
			"AA******AA********A.....",
			"..CC****CCCC******CC....",
			"....GG**GGGGGG****GGG...",
			"......TTTTTTTTTT**TTTT.."};
		String[] taxons = new String[]{
			"srpA","srpC","srpG","srpT"
		};
		
		for (int i = 0; i < srpMap.getSrpCount(); i++) {
			Sequence s = new Sequence(new Taxon(taxons[i]), seqs[i]);
			ShortRead srp = new ShortRead(s);
			assertEquals(srp.getStart(), srpMap.getSrpStart(i));
			assertEquals(srp.getEnd(), srpMap.getSrpEnd(i));
			assertEquals(srp.getFragmentSrp(), srpMap.getSrpFragment(i));
			assertEquals(srp.getFullSrp(), srpMap.getSrpFull(i));
			assertEquals(srp.getName(), srpMap.getSrpName(i));
			assertEquals(srp.getLength(), srpMap.getSrpLength(i));
		}
		
		
	}
	

	@Test
	public void testMapNameToID() {
		String[] taxons = new String[]{
				"srpA","srpC","srpG","srpT"
		};
		for (int i = 0; i < taxons.length; i++) {
			Integer actual= srpMap.mapNameToID(taxons[i]);
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
			assertEquals(expected, srpMap.getMapToSrp(2*i));
			assertEquals(expected, srpMap.getMapToSrp(2*i+1));
		}
		for (int i = 0; i < 4; i++) {
			assertEquals(expected, srpMap.getMapToSrp(18+i));
			expected.remove(0);
		}
		assertEquals(expected, srpMap.getMapToSrp(22));
	}
}
