package test.srp.operator.haplotypes;


import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.evolution.OperationRecord;
import srp.evolution.OperationType;
import srp.haplotypes.Haplotype;
import srp.haplotypes.HaplotypeModel;
import srp.operator.haplotypes.AbstractHaplotypeOperator;
import srp.operator.haplotypes.BaseSingleOperator;
import test.TestUtils;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;

public class BaseSingleOperatorTest {

	private static DataType datatype = Nucleotides.INSTANCE;
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

		AbstractHaplotypeOperator operator = new BaseSingleOperator(haplotypeModel);
    	assertEquals(operator.getOperatorName(), "BaseSingleOperator");
    	assertEquals(operator.getPerformanceSuggestion(), "");
    	assertEquals(operator.getOperationType(), OperationType.SINGLE);
	}

	@Test
	public void testOperator() throws Exception {
//		String[] seqs = new String[]{
//				"AACCGGTT",
//				"...GCTAT",
//				"ACCGT...",
//				};

		HaplotypeModel haplotypeModel = new HaplotypeModel(5, 10);//AlignmentUtils.createAlignment(seqs));
		SimpleMCMCOperator op = new BaseSingleOperator(haplotypeModel);

		char[][] storedChars = new char[haplotypeModel.getHaplotypeCount()][haplotypeModel.getHaplotypeLength()];
		for (int i = 0; i < storedChars.length; i++) {
			Haplotype haplotype = haplotypeModel.getHaplotype(i);
			haplotype.getChars(0, haplotype.getLength(), storedChars[i], 0);
		}
		
		int swapHapCount[] = new int[haplotypeModel.getHaplotypeCount()];
				
		int ite = 100000;
		for (int o = 0; o < ite; o++) {
			try {
				op.doOperation();
				
				OperationRecord opRecord = haplotypeModel.getOperationRecord();
				int hapIndex = opRecord.getSpectrumIndex();
				int siteIndex = opRecord.getSingleIndex();

				char newChar = haplotypeModel.getHaplotype(hapIndex).getChar(siteIndex);
				char oldChar = storedChars[hapIndex][siteIndex];

				swapHapCount[hapIndex]++;
				assertNotEquals(newChar, oldChar);
				int newState = datatype.getState(newChar);
				assertTrue(newState<4);
				
				storedChars[hapIndex][siteIndex] = newChar;
				
			} catch (OperatorFailedException e) {
//				e.printStackTrace();
			}	
		}	
		
		double expectedMean = 1.0/swapHapCount.length;
		for (int i = 0; i < swapHapCount.length; i++) {
			TestUtils.assertExpectationRange(swapHapCount[i]/(double) ite, expectedMean, 0.005);  
		}
	}

}
