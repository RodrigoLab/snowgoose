package test.alignment;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.nio.file.Files;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import alignment.AlignmentMapping;
import alignment.ShortRead;
import core.DataImporter;
import dr.evolution.alignment.Alignment;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

public class AlignmentMappingTest {

	// private static Alignment shortReads;
	private static DataImporter dataImporter;
	private static String refFileName;
	private static String srpFileName;
	private Alignment shortReads;
//	private Sequence refSeq;

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		
		String dir = System.getProperty("user.dir")+File.separatorChar+"data"+File.separator;
		System.out.println(dir);
//		/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/";
//		refFileName = "1110_10.ref";
		srpFileName = "align_test_4.fasta";
		dataImporter = new DataImporter(dir);

	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {

		shortReads = dataImporter.importAlignment(srpFileName);
//		refSeq = dataImporter.importRefSeq(refFileName);

	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testAlignmentMapping() {

		System.out.println("===");
		System.out.println(shortReads.getPatternCount());
		System.out.println(shortReads.getSiteCount());
		System.out.println(shortReads.getStateCount());
		System.out.println(shortReads.getTaxonCount());
		// System.out.println(shortReads.get
		for (int i = 0; i < shortReads.getSequenceCount(); i++) {
			Sequence s = shortReads.getSequence(i);
			Taxon t = s.getTaxon();

		}

		AlignmentMapping amap = new AlignmentMapping(shortReads);
		
		ShortRead srp = amap.getShortRead(0);
		String expected = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA...........................................................................";
		assertEquals(expected, srp.getFullSrp());
		expected = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
		assertEquals(expected, srp.getFragmentSrp());
		
	}

	@Test
	public void testAddSequence() {
		System.out.println(shortReads.getId());
	}

}
