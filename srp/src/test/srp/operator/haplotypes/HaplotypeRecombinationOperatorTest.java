package test.srp.operator.haplotypes;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.evolution.OperationType;
import srp.evolution.haplotypes.HaplotypeModel;
import srp.operator.haplotypes.AbstractHaplotypeOperator;
import srp.operator.haplotypes.HaplotypeRecombinationOperator;

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
		
		HaplotypeModel haplotypeModel = new HaplotypeModel(5, 100);

		AbstractHaplotypeOperator operator = new HaplotypeRecombinationOperator(haplotypeModel, 0);
    	assertEquals(operator.getOperatorName(), "HaplotypeRecombinationOperator");
    	
    	assertEquals(operator.getOperationType(), OperationType.RECOMBINATION);
	}
	
	
	@Test
	public void testOp() throws Exception {
		HaplotypeModel haplotypeModel = new HaplotypeModel(2, 100);
		AbstractHaplotypeOperator operator = new HaplotypeRecombinationOperator(haplotypeModel, 0);
		
		for (int i = 0; i < 1000; i++) {
			
			String expected0 = haplotypeModel.getAlignedSequenceString(0);
			String expected1 = haplotypeModel.getAlignedSequenceString(1);
			
			operator.doOperation();
			
			int[] rPos = haplotypeModel.getOperationRecord().getRecombinationPositionIndex();
			String hap0 = haplotypeModel.getAlignedSequenceString(0);
			String hap1 = haplotypeModel.getAlignedSequenceString(1);
			
			assertNotEquals(expected0, hap0);
			assertNotEquals(expected1, hap1);
			
			assertEquals(expected0.substring(0, rPos[0]), hap0.substring(0, rPos[0]));
			assertEquals(expected1.substring(0, rPos[0]), hap1.substring(0, rPos[0]));
			
			assertEquals(expected0.substring(rPos[0]), hap1.substring(rPos[0]));
			assertEquals(expected1.substring(rPos[0]), hap0.substring(rPos[0]));
			
			
		}
		
		
		
    	
	}
	
}
