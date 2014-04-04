package test.srp.operator.haplotypes;


import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.apache.commons.lang3.StringUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.evolution.haplotypes.old.OldHaplotypeModel;
import srp.evolution.haplotypes.old.OldHaplotypeModelUtils;
import srp.evolution.shortreads.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.operator.haplotypes.HaplotypeRecombinationOperator;
import dr.inference.operators.SimpleMCMCOperator;


public class HaplotypeRecombinationOperatorTest {

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
	public void testGetOperatorName() {
		String[] seqs = new String[]{
				"AAAAAAAAATGTGTTTT....",
				".....CCCCCCCCCCCCCCCCCCCTTTTCCCC....",
				"..........GGGGGGGGGGGGGGCGCGTATAGGGG",
				"...............TTTTTTTTTACACTATA....",
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		OldHaplotypeModel haplotypeModel = new OldHaplotypeModel(aMap, 3);

		SimpleMCMCOperator operator = new HaplotypeRecombinationOperator(haplotypeModel, 0);
    	assertEquals(operator.getOperatorName(), "HaplotypeRecombinationOperator");
    	assertEquals(operator.getPerformanceSuggestion(), "");
    	
	}

	
	@Test
	public void testDoOperation() throws Exception {
		
		String[] srp = new String[]{
				".....TTTTT",
				"TTTTT....."
				};
		
		String[] haps = new String[]{
				"AAAAAAAAAA",
				"CCCCCCCCCC"
				};
		
		for (int i = 1; i < 1000; i++) {
			OldHaplotypeModel haplotypeModel = OldHaplotypeModelUtils.createHaplotypeModel(srp, haps);
			HaplotypeRecombinationOperator op = new HaplotypeRecombinationOperator(haplotypeModel, i);
			op.doOperation();
			String s0 = haplotypeModel.getHaplotypeString(0);
			String s1 = haplotypeModel.getHaplotypeString(1);
			assertTrue( StringUtils.countMatches(s0, "A")>0 && StringUtils.countMatches(s0, "A")<10);
			assertTrue( StringUtils.countMatches(s0, "C")>0 && StringUtils.countMatches(s0, "C")<10);
			assertTrue( StringUtils.countMatches(s1, "A")>0 && StringUtils.countMatches(s1, "A")<10);
			assertTrue( StringUtils.countMatches(s1, "C")>0 && StringUtils.countMatches(s1, "C")<10);
			assertEquals(10, StringUtils.getLevenshteinDistance(s0, s1));

			haplotypeModel.reject();
			assertEquals(haps[0], haplotypeModel.getHaplotypeString(0));
			assertEquals(haps[1], haplotypeModel.getHaplotypeString(1));

		}
	}
		
	
	

}
