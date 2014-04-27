package test.srp.operator.haplotypes;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.evolution.OperationType;
import srp.haplotypes.HaplotypeModel;
import srp.operator.haplotypes.AbstractHaplotypeOperator;
import srp.operator.haplotypes.HaplotypeSwapSectionOperator;
import dr.math.MathUtils;

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
		
		HaplotypeModel haplotypeModel = new HaplotypeModel(5, 100);

		AbstractHaplotypeOperator operator = new HaplotypeSwapSectionOperator(haplotypeModel, 0, null);
    	assertEquals(operator.getOperatorName(), "HaplotypeSwapSectionOperator");
    	
    	assertEquals(operator.getOperationType(), OperationType.RECOMBINATION);
	}
	
	
	@Test
	public void testOp() throws Exception {
		HaplotypeModel haplotypeModel = new HaplotypeModel(2, 100);
		AbstractHaplotypeOperator operator = new HaplotypeSwapSectionOperator(haplotypeModel, 10, null);
		
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
			
			assertEquals(expected0.substring(rPos[0], rPos[1]), hap1.substring(rPos[0], rPos[1]));
			assertEquals(expected1.substring(rPos[0], rPos[1]), hap0.substring(rPos[0], rPos[1]));
			
			assertEquals(expected0.substring(rPos[1]), hap0.substring(rPos[1]));
			assertEquals(expected1.substring(rPos[1]), hap1.substring(rPos[1]));
			
			assertEquals(10, rPos[1]-rPos[0]);
			
		}
	}

	@Test
	public void testOpMultiH() throws Exception {
		
		int hapCount = 10;
		String[] expecteds = new String[hapCount];
		HaplotypeModel haplotypeModel = new HaplotypeModel(10, 200);
		AbstractHaplotypeOperator operator = new HaplotypeSwapSectionOperator(haplotypeModel, 10, null);
		int bases = 0;
		for (int i = 0; i < 1000; i++) {
			do{
				bases = MathUtils.nextInt(50);
			}while(bases==0);

			operator = new HaplotypeSwapSectionOperator(haplotypeModel, bases, null);

			for (int h = 0; h < hapCount; h++) {
				expecteds[h] =  new String(haplotypeModel.getAlignedSequenceString(h));
			}
			
			operator.doOperation();
			
			int[] rPos = haplotypeModel.getOperationRecord().getRecombinationPositionIndex();
			int[] hIndex = haplotypeModel.getOperationRecord().getRecombinationSpectrumIndex();
			
			int count = 0;
			for (int h = 0; h < hapCount; h++) {
				if( h!= hIndex[0] && h!=hIndex[1] ){
					assertEquals(expecteds[h], haplotypeModel.getAlignedSequenceString(h));
					count++;
				}
			}
			assertEquals(hapCount-2, count);
			
			String hap0 = haplotypeModel.getAlignedSequenceString(hIndex[0]);
			String hap1 = haplotypeModel.getAlignedSequenceString(hIndex[1]);
		
			assertEquals(expecteds[hIndex[0]].substring(0, rPos[0]), hap0.substring(0, rPos[0]));
			assertEquals(expecteds[hIndex[1]].substring(0, rPos[0]), hap1.substring(0, rPos[0]));
			
			assertEquals(expecteds[hIndex[0]].substring(rPos[0], rPos[1]), hap1.substring(rPos[0], rPos[1]));
			assertEquals(expecteds[hIndex[1]].substring(rPos[0], rPos[1]), hap0.substring(rPos[0], rPos[1]));
			
			assertEquals(expecteds[hIndex[0]].substring(rPos[1]), hap0.substring(rPos[1]));
			assertEquals(expecteds[hIndex[1]].substring(rPos[1]), hap1.substring(rPos[1]));
			
			assertEquals(bases, rPos[1]-rPos[0]);
			
		}
	}

}
