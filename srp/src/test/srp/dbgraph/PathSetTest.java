package test.srp.dbgraph;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import com.google.common.primitives.Ints;

import srp.dbgraph.DeBruijnImporter;
import srp.dbgraph.Path;
import srp.dbgraph.PathSet;

public class PathSetTest {

	private String dataDir;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/DBGraphData/";
	}

	@After
	public void tearDown() throws Exception {
	}

	
	
	@Test
	public void testImport() throws Exception {
		String fileName = "H7_50.60.pathsc.txt";
		PathSet pathSets = DeBruijnImporter.importPathSets(dataDir, fileName, true);
		System.out.println(pathSets.getPathCount());
		Path path = pathSets.getPath(1);
		ArrayList<Integer> nodeList = path.getNodeList();
		int[] array = Ints.toArray(nodeList);
		
		
		int[] expecteds = new int[]{6, 7, 14, 21, 19, 27, 34, 45, 46, 53, 47, 56, 62, 63, 71, 76, 85, 92, 103, 109, 118, 124, 131};
		assertArrayEquals(expecteds, array);
	}
}
