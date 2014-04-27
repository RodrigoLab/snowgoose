package test.srp.dbgraph;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.dbgraph.CompatibleNode;
import srp.dbgraph.CompatibleSets;

public class CompatibleSetsTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
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
	public void testConstructor() throws Exception {
		
		CompatibleSets cSets = new CompatibleSets();
		CompatibleNode cNode = new CompatibleNode(
				"4:30 4:0:4903 18:31:554 30:32:20056 17:61:3745 27:66:28245 41:126:106 55:128:55372 63:405:3356 72:429:218 82:431:3350");
		cSets.addCompatibleNode(cNode);

		cNode = new CompatibleNode("11:31 11:0:4060 18:31:462 29:32:97682 ");
		cSets.addCompatibleNode(cNode);

		cNode = new CompatibleNode(
				"18:31 30:32:1038 29:32:6645 17:61:220 27:66:1831 41:126:4 55:128:3302 63:405:161 72:429:11 82:431:325");
		cSets.addCompatibleNode(cNode);

		cNode = new CompatibleNode(
				"12:39 12:0:11591 22:40:10538 38:49:72912 46:229:8899 61:263:281 19:264:1103 31:267:22140 43:308:74195 58:513:55 67:514:745");
		cSets.addCompatibleNode(cNode);

		assertEquals(18, cSets.getMaxNodeIndex());
		
		cSets.preprocess();
		assertEquals(4, cSets.getNodeCount());
		int[] expecteds = new int[]{4, 11, 18, 12};
		int[] expectedsDepth = new int[]{30, 31, 31, 39};
		int[] expectedsCNodeLength = new int[]{10, 3, 9, 10};
		for (int i = 0; i < cSets.getNodeCount(); i++) {
			assertEquals(expecteds[i], cSets.getCompatibleNode(i).getNodeIndex() );
			assertEquals(expectedsDepth[i], cSets.getCompatibleNode(i).getNodeDepth() );
			assertEquals(expectedsCNodeLength[i], cSets.getCompatibleNode(i).getCNodeArray().length );
			
		}
	}
}
