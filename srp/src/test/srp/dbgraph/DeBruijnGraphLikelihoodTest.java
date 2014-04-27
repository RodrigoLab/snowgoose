package test.srp.dbgraph;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.dbgraph.CompatibleSets;
import srp.dbgraph.DeBruijnGraph;
import srp.dbgraph.DeBruijnGraphLikelihood;
import srp.dbgraph.DeBruijnImporter;

public class DeBruijnGraphLikelihoodTest {

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
	public void testBasic() throws Exception {
		DeBruijnImporter dbi = new DeBruijnImporter(dataDir);
		DeBruijnGraph dbg = dbi.importDeBruijnGraph("N10_cond.graph");
		CompatibleSets cSets = dbi.importCompatibleSet("N10_comp.txt");
		
		DeBruijnGraphLikelihood dbgLikelihood = new DeBruijnGraphLikelihood(dbg, cSets);
	}
}
