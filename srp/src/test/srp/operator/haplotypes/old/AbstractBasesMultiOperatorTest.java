package test.srp.operator.haplotypes.old;


import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.evolution.haplotypes.old.OldHaplotypeModel;
import srp.evolution.shortreads.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.operator.haplotypes.old.AbstractBasesMultiOperator;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

public class AbstractBasesMultiOperatorTest {

	private OldHaplotypeModel haplotypeModel;

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		String[] seqs = new String[]{
				"AAAAAAAAATGTGTTTT....",
				".....CCCCCCCCCCCCCCCCCCCTTTTCCCC....",
				"..........GGGGGGGGGGGGGGCGCGTATAGGGG",
				"...............TTTTTTTTTACACTATA....",
				"CCCCCTTTTTAAAAAGGGGGTCGATGCAGTAGCTAG"
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		haplotypeModel = new OldHaplotypeModel(aMap, 3);

	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testGetCoercableParameter() throws Exception {

		int nBases = 10;
		CoercableMCMCOperator operator = new MockAbstractSwapBaseseOperator(haplotypeModel, nBases, null);

    	
    	assertEquals(operator.getRawParameter(), nBases, 0);
    	assertEquals(operator.getCoercableParameter(), Math.log(nBases-1), 1e-10); 
    	
		double autoOptimise = 20;
		operator.setCoercableParameter(autoOptimise);
		assertEquals(operator.getCoercableParameter(), autoOptimise, 0);
		assertEquals(operator.getRawParameter(), haplotypeModel.getHaplotypeLength(), 1e-10);

		autoOptimise = 2;
		operator.setCoercableParameter(autoOptimise);
		assertEquals(operator.getCoercableParameter(), autoOptimise, 0);
		assertEquals(operator.getRawParameter(), (int) Math.exp(autoOptimise)+1, 0);

		autoOptimise = -10;
		operator.setCoercableParameter(autoOptimise);
		assertEquals(operator.getCoercableParameter(), autoOptimise, 0);
		assertEquals(operator.getRawParameter(), 1, 0);
    	
    	
	}
	
	@Test
	public void testGetCoercableParameter2() throws Exception {

		String[] seqs = new String[]{
				".......................................................................................................................................................................................................A",				
				"CCCCCTTTTTAAAAAGGGGG....."
				};
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(seqs);
		haplotypeModel = new OldHaplotypeModel(aMap, 3);

		int nBases = 10;
		CoercableMCMCOperator operator = new MockAbstractSwapBaseseOperator(haplotypeModel, nBases, null);
		final int scaleFactor = (int) (haplotypeModel.getHaplotypeLength()*0.01);
    	
    	assertEquals(operator.getRawParameter(), nBases, 0);
    	assertEquals(operator.getCoercableParameter(), Math.log(nBases-1)/scaleFactor, 1e-10); 
    	
		double autoOptimise = 20;
		operator.setCoercableParameter(autoOptimise);
		assertEquals(operator.getCoercableParameter(), autoOptimise, 0);
		assertEquals(operator.getRawParameter(), (int) Math.exp(autoOptimise*scaleFactor)+1, 1e-10);

		autoOptimise = 2;
		operator.setCoercableParameter(autoOptimise);
		assertEquals(operator.getCoercableParameter(), autoOptimise, 0);
		assertEquals(operator.getRawParameter(), (int) Math.exp(autoOptimise*scaleFactor)+1, 0);

		autoOptimise = -10;
		operator.setCoercableParameter(autoOptimise);
		assertEquals(operator.getCoercableParameter(), autoOptimise, 0);
		assertEquals(operator.getRawParameter(), 1, 0);
    	
    	
		
	}
	
	@Test
	public void testResetAllPosChars() throws Exception {
	
		int nBases = 10;
		MockAbstractSwapBaseseOperator operator = new MockAbstractSwapBaseseOperator(haplotypeModel, nBases, null);
		final int scaleFactor = (int) (haplotypeModel.getHaplotypeLength()*0.01);
    	
		int[][] allPosChars = operator.testGetAllPosChars();
		int[] expecteds = new int[haplotypeModel.getHaplotypeLength()];
		assertArrayEquals(expecteds , allPosChars[0]);
		assertArrayEquals(expecteds, allPosChars[1]);
		
		Arrays.fill(allPosChars[0], 1);
		Arrays.fill(allPosChars[1], 2);
		operator.testSetAllPosChars(allPosChars);

		Arrays.fill(expecteds, 1);
		assertArrayEquals(expecteds , allPosChars[0]);
		Arrays.fill(expecteds, 2);
		assertArrayEquals(expecteds, allPosChars[1]);
		
		operator.testResetAllPosChars();
		Arrays.fill(expecteds, -1);
		assertArrayEquals(expecteds, allPosChars[0]);
		assertArrayEquals(expecteds, allPosChars[1]);
		
	}
	
	private class MockAbstractSwapBaseseOperator extends AbstractBasesMultiOperator {

		public MockAbstractSwapBaseseOperator(OldHaplotypeModel haplotypeModel,
				int swapLength, CoercionMode mode) {
			super(haplotypeModel, swapLength, mode);

		}

		public void testResetAllPosChars(){
			resetAllPosChars();
		}
		
		public int[][] testGetAllPosChars(){
			return allPosChars;
		}
		
		public void testSetAllPosChars(int[][] intArrays){
			allPosChars = intArrays;
		}

		@Override
		public String getPerformanceSuggestion() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public String getOperatorName() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public double doOperation() throws OperatorFailedException {
			// TODO Auto-generated method stub
			return 0;
		}
	}
}
