package test.srp.dbgraph;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.dbgraph.DeBruijnImporter;

public class DeBruijnImporterTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}


	private String dataDir;

	@Before
	public void setUp() throws Exception {
		dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/DBGraphData";
	}

	@After
	public void tearDown() throws Exception {
	}

	
	@Test
	public void testImportDeBruijnGraph() throws Exception {
		DeBruijnImporter dbi = new DeBruijnImporter(dataDir);
		dbi.importDeBruijnGraph("H7_56-Paired.bfast.60.cond.graph");
	}
	

	@Test
	public void testImportCompatibleSet() throws Exception {
		DeBruijnImporter dbi = new DeBruijnImporter(dataDir);
		dbi.importCompatibleSet("H7_56-Paired.bfast.60.comp.txt");
	}
}
