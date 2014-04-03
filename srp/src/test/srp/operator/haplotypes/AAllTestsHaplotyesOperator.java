package test.srp.operator.haplotypes;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ 
	AbstractBasesMultiOperatorTest.class,
	
	ColumnOperatorTest.class,
	
	HaplotypeRecombinationOperatorTest.class,
	HaplotypeSwapSectionOperatorTest.class,
	
	BasesMultiEmpiricalOperatorTest.class,
	BasesMultiOperatorTest.class,
	BasesMultiUniformOperatorTest.class,
	
	BaseSingleFrequencyOperatorTest.class,
	BaseSingleOperatorTest.class, 
	BaseSingleUniformOperatorTest.class,


})
public class AAllTestsHaplotyesOperator {

}
