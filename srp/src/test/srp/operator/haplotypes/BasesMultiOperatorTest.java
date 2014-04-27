package test.srp.operator.haplotypes;

import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.evolution.OperationType;
import srp.haplotypes.HaplotypeModel;
import srp.operator.haplotypes.AbstractHaplotypeOperator;
import srp.operator.haplotypes.BasesMultiOperator;

public class BasesMultiOperatorTest {

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

		AbstractHaplotypeOperator operator = new BasesMultiOperator(haplotypeModel, 5, null);
    	assertEquals(operator.getOperatorName(), "BasesMultiOperator");
    	assertEquals(operator.getOperationType(), OperationType.MULTI);
    	assertEquals(operator.getRawParameter(), 5, 0);
	}
	
}
