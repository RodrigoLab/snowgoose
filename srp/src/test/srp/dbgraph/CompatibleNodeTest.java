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
		String input = "4:30 4:0:4903 18:31:554 30:32:20056 17:61:3745 27:66:28245 41:126:106 55:128:55372 63:405:3356 72:429:218 82:431:3350";
		
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
		assertEquals(4, compNode.getNodeIndex());
		assertEquals(30, compNode.getDepth());
		compNode.getCNodeArray();
		
	}
}
