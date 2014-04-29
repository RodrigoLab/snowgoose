package test.srp.dbgraph;

import static org.junit.Assert.*;

import java.util.Arrays;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import com.google.common.primitives.Ints;

import srp.dbgraph.CompatibleSets;
import srp.dbgraph.DeBruijnGraph;
import srp.dbgraph.DeBruijnGraphLikelihood;
import srp.dbgraph.DeBruijnImporter;
import srp.dbgraph.Path;
import srp.dbgraph.PathSet;

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
		double[][] lengthDiff = dbgLikelihood.getLengthDiff();
		
		double[][] expecteds = new double[15][15];
		expecteds[0][10] = 8 - 8 + 40 + 0 + 10;
		expecteds[1][9] = 10 - 10 + 40 + 1 + 9;
		expecteds[2][10] = 15 - 15 + 40 + 2 + 10;
		expecteds[3][11] = 6 - 60 + 40 + 3+ 11;
		expecteds[4][12] = 3-30 + 40 + 4 + 12;
		expecteds[4][11] = 6 - 30 + 40 + 4 + 11;
		expecteds[6][13] = 26 - 26 + 40 + 6 + 13;
		expecteds[7][14] = 16 - 16 +40 + 7 + 14;
		expecteds[9][10] = 10 - 18 + 40 + 9 + 10;
		expecteds[11][12] = 31 - 31 + 40 + 11 + 12;
		expecteds[12][13] = 26 - 39 + 40  + 12 + 13;
		expecteds[13][8] = 83- 40 + 40 + 13 + 8;
		expecteds[13][10] = 15 - 40 + 40 + 13 + 10; 
		expecteds[13][14] = 16 - 40 + 40 + 13 + 14;
		expecteds[14][0] = 0 - 85 + 40 + 14 + 0;
		expecteds[14][8] = 83 - 85 + 40 + 14 + 8;
		expecteds[14][10] = 15 - 85 + 40 + 14 + 10;
		for (int i = 0; i < expecteds.length; i++) {
			assertArrayEquals("node i error? "+i, expecteds[i], lengthDiff[i], 0);
		}
		int expectedCount = 776;
		assertEquals(expectedCount, dbgLikelihood.getTotalCount());
		

//		0:8 0:0:128 10:8:14
//		1:10 1:0:25 9:10:47
//		2:15 2:0:39 10:15:37
//		3:60 3:0:50 11:6:133
//		4:30 4:0:49 12:3:55 11:6:37
//		5:11 5:0:15
//		6:26 6:0:11 13:26:26
//		7:16 7:0:49 14:16:38
//		8:84 8:0:122
//		9:18 9:0:21 10:10:44
//		10:19 10:0:157
//		11:31 11:0:40 12:31:46 
//		12:39 12:0:11 13:26:11
//		13:40 13:0:10 8:83:31 10:15:83 14:16:74
//		14:85 0:0:9 8:83:9 10:15:82 14:16:52
	}
	
	
	@Test
	public void testLikelihood() throws Exception {
		
		DeBruijnImporter dbi = new DeBruijnImporter(dataDir);
		DeBruijnGraph dbGraph = dbi.importDeBruijnGraph("H7_50-Paired.bfast.60.cond.graph");
		CompatibleSets compSets = dbi.importCompatibleSet("H7_50-Paired.bfast.60.comp.txt");
		PathSet allPaths = dbi.importPathSets("H7_50.60.pathsc.txt", true);
//				H7_50.60.like.txt    H7_50-Paired.bfast.60.comp.txt
//				H7_50.60.pathsc.txt  H7_50-Paired.bfast.60.cond.graph
		for (int i = 1; i < allPaths.getPathCount(); i++) {
			
			Path path = allPaths.getPath(i);
//			System.out.println(path.getNodeList());
		}
		DeBruijnGraphLikelihood dbgLikelihood= new DeBruijnGraphLikelihood(dbGraph, compSets);
//		double ln = dbgLikelihood.calculateLikelihood(allPaths);
//		System.out.println(ln);
		
		double expected;
		double ln;
		int[] nodeList = new int[]{63,26,119,30,27,25,28,75,14,92,24,31,78,121,22,106,13,23,29,6,3,36,51,34,37,62};
		PathSet subPaths = new PathSet(false);
		for (int i : nodeList) {
			Path path = allPaths.getPath(i);
//			System.out.println(i +"\t"+ path.getNodeList());
			subPaths.addPath(path);
		}
//		System.out.println(subPaths.);
//		63 26 119 30 27 25 28 75 14 92 24 31 78 121 22 106 13 23 29 6 3 36 51 34 37 62 -728496975.229625
//		63 7 119 30 27 25 75 14 92 24 31 91 78 121 22 106 13 23 29 6 3 36 51 34 37 62 -731239363.747458
//		63 26 119 30 27 25 28 75 24 22 106 13 29 6 36 37 -727686676.708421
//		63 119 30 27 25 75 24 22 106 13 29 6 36 37 -725722248.160968
//		expected = -728496975.229625;
//		ln = dbgLikelihood.calculateLikelihood(subPaths);
//		System.out.println(ln);
//		assertEquals(expected, ln, 1e-8);
		
		
		nodeList = new int[]{106,119,13,22,24,25,27,29,30,36,37,6,63,75};
//		Arrays.sort(nodeList);
//		ArrayUtils.reverse(nodeList);
		subPaths = new PathSet(false);
		for (int i : nodeList) {
			Path path = allPaths.getPath(i);
//			System.out.println(path.getNodeCount());
			subPaths.addPath(path);
		}
		
		expected = -725722248.160968;
		ln = dbgLikelihood.calculateLikelihood(subPaths);
		System.out.println(ln);
		assertEquals(expected, ln, 1e-8);
		
			
		double[][] computeDHashTable = dbgLikelihood.computeDHashTable(subPaths);
		for (int i = 0; i < computeDHashTable.length; i++) {
//			System.out.println(Arrays.toString(computeDHashTable[i]));
//			double sum = StatUtils.sum(computeDHashTable[i]);
//			System.out.println(i +"\t"+ sum);
		}
		
	}
}
