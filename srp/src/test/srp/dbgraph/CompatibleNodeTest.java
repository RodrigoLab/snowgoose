package test.srp.dbgraph;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.dbgraph.CompatibleNode;

import com.google.common.base.Splitter;

public class CompatibleNodeTest {

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
		String input = "18:31 30:32:1038 29:32:6645 17:61:220 27:66:1831 41:126:4 55:128:3302 63:405:161 72:429:11 82:431:325";
		
		String regex = "\\d+\\s\\d+";
		
		Iterable<String> split = Splitter.on(" ").split(input.trim());
		CompatibleNode compNode = null;//new CompatibleSet();
		for (String string : split) {
			if(compNode == null){
				Iterable<String> split2 = Splitter.on(":").split(string);
				compNode = new CompatibleNode(split2);
			}
			else{
				compNode.addCompatibleNode(string);
			}
		}
//			compatibleSets.addCompatibleNode(compNode);
//			
//		CompatibleNode cNode = new CompatibleNode();
		compNode.proprocess();
		assertEquals(18, compNode.getNodeIndex());
		assertEquals(31, compNode.getNodeDepth());
		int[] expecteds = new int[]{30,29,17,27,41,55,63,72,82};
		assertArrayEquals(expecteds, compNode.getCNodeArray() );
		
		int[][] expecteds2D = new int[][]{
				{30,32,1038},
				{29,32,6645},
				{17,61,220},
				{27,66,1831},
				{41,126,4},
				{55,128,3302},
				{63,405,161},
				{72,429,11},
				{82,431,325}};
		
		for (int i = 0; i < expecteds2D.length; i++) {
			assertEquals(expecteds2D[i][1] ,compNode.getCNodeDepth(expecteds2D[i][0]) );
			assertEquals(expecteds2D[i][2] ,compNode.getCNodeCount(expecteds2D[i][0]) );
		}
		
	}
	
	@Test
	public void testConstructorString() throws Exception {
		String input = "4:30 4:0:4903 18:31:554 30:32:20056 17:61:3745 27:66:28245 41:126:106 55:128:55372 63:405:3356 72:429:218 82:431:3350";
		CompatibleNode compNode = new CompatibleNode(input);
		
		compNode.proprocess();
		assertEquals(4, compNode.getNodeIndex());
		assertEquals(30, compNode.getNodeDepth());
		int[] expecteds = new int[]{4, 18, 30, 17, 27, 41, 55, 63, 72, 82};
		assertArrayEquals(expecteds, compNode.getCNodeArray() );
		
		int[][] expecteds2D = new int[][]{
				{4,0,4903},
				{18,31,554},
				{30,32,20056},
				{17,61,3745},
				{27,66,28245},
				{41,126,106},
				{55,128,55372},
				{63,405,3356},
				{72,429,218},
				{82,431,3350}};
		
		for (int i = 0; i < expecteds2D.length; i++) {
			assertEquals(expecteds2D[i][1] ,compNode.getCNodeDepth(expecteds2D[i][0]) );
			assertEquals(expecteds2D[i][2] ,compNode.getCNodeCount(expecteds2D[i][0]) );
		}
	}
}
