package test.srp.haplotypes.operator;


import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.operator.UniformSwapBaseOperator;
import dr.evolution.alignment.SimpleAlignment;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;

public class UniformSwapBaseOperatorTest {

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
		
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, 3);

		MCMCOperator operator = new UniformSwapBaseOperator(haplotypeModel, 0);
    	assertEquals(operator.getOperatorName(), "UniformSwapBaseOperator");
	}

	@Test
	public void testDoOperation() throws OperatorFailedException {
		String[] seqs = new String[]{
				"AAA.",
				"AAA.",
				"AAA.",
				"AAA.",
				"AAA.",
				"AAA.",
				"AAA.",
				"A*A.",
				"AAC.",
				"AAG.",
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		
		String[] haps = new String[]{
				"TTTT"
				};
		SimpleAlignment hapAlignment = AlignmentUtils.createAlignment(haps);
		HaplotypeModel haplotypeModel = new HaplotypeModel(aMap, hapAlignment);

		SimpleMCMCOperator operator = new UniformSwapBaseOperator(haplotypeModel, 0);
    	
    	
    	for (int i = 0; i < 100; i++) {
    		operator.doOperation();
			String newHap = haplotypeModel.getHaplotypeString(0);
			assertNotEquals(haps[0], newHap);
			int pos = haplotypeModel.getSwapInfo().getSwapInfoSWAPBASE()[1];
			assertTrue( newHap.charAt(pos)<'T');
//			haps[0] = newHap;
			haplotypeModel.reject();
			newHap = haplotypeModel.getHaplotypeString(0);
			assertEquals(haps[0], newHap);

		}

	}
}
