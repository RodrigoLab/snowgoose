package test.srp.operator.haplotypes;


import static org.junit.Assert.assertEquals;

import org.apache.commons.lang3.StringUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.old.OldHaplotypeModel;
import srp.haplotypes.old.OldHaplotypeModelUtils;
import srp.operator.haplotypes.HaplotypeSwapSectionOperator;
import srp.shortreads.AlignmentMapping;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.SimpleMCMCOperator;

public class HaplotypeSwapSectionOperatorTest {

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

		int nBases = 10;
		CoercableMCMCOperator operator = new HaplotypeSwapSectionOperator(haplotypeModel, nBases, null);
    	assertEquals(operator.getOperatorName(), "HaplotypeSwapSectionOperator");
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


		for (int i = 1; i < 10; i++) {
			OldHaplotypeModel haplotypeModel = OldHaplotypeModelUtils.createHaplotypeModel(srp, haps);
			SimpleMCMCOperator op = new HaplotypeSwapSectionOperator(haplotypeModel, i, CoercionMode.COERCION_OFF);

			for (int j = 0; j < 100; j++) {
				op.doOperation();
				String s0 = haplotypeModel.getHaplotypeString(0);
				String s1 = haplotypeModel.getHaplotypeString(1);
				assertEquals(10-i, StringUtils.countMatches(s0, "A"));
				assertEquals(i, StringUtils.countMatches(s0, "C"));
				assertEquals(i, StringUtils.countMatches(s1, "A"));
				assertEquals(10-i, StringUtils.countMatches(s1, "C"));
				haplotypeModel.reject();
				assertEquals(haps[0], haplotypeModel.getHaplotypeString(0));
				assertEquals(haps[1], haplotypeModel.getHaplotypeString(1));
			}
		}
		
	}
	
}