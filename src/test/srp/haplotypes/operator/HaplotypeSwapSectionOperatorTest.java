package test.srp.haplotypes.operator;


import static org.junit.Assert.assertEquals;

import org.apache.commons.lang3.StringUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.operator.HaplotypeSwapSectionOperator;

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
			HaplotypeModel haplotypeModel = AlignmentUtils.createHaplotypeModel(srp, haps);
			HaplotypeSwapSectionOperator op = new HaplotypeSwapSectionOperator(haplotypeModel, i, null);

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