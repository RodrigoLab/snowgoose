package test.alignment;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import core.DataImporter;
import dr.evolution.alignment.Alignment;

import alignment.AlignmentMapping;
import alignment.AlignmentMatrix;

public class AlignmentMatrixTest {

	private static AlignmentMapping amap;

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		
		String dir = "/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/";
		String refFileName = "1110_10.ref";
		String srpFileName = "1110_10_align_test.fasta";
		DataImporter dataImporter = new DataImporter(dir);

		Alignment shortReads = dataImporter.importAlignment(srpFileName);
		amap = new AlignmentMapping(shortReads);
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
	public void testAlignmentMatrix() {
	
		AlignmentMatrix alignmentMatrix = new AlignmentMatrix(amap, 5);
		alignmentMatrix.testGetSeq();

		alignmentMatrix.testMultipleSeq();
		alignmentMatrix.toAlignment();
	}
	@Test
	public void testRandomSeq() {

	}

	@Test
	public void testSwapSrp() {
		fail("Not yet implemented");
	}

}
