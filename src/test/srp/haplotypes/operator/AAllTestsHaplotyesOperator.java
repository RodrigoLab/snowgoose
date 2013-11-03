package test.srp.haplotypes.operator;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ 
	AbstractMultiBasesOperatorTest.class,
	
	ColumnOperatorTest.class,
	
	HaplotypeRecombinationOperatorTest.class,
	HaplotypeSwapSectionOperatorTest.class,
	
	MultiBasesEmpiricalOperatorTest.class,
	MultiBasesOperatorTest.class,
	MultiBasesUniformOperatorTest.class,
	
	SingleBaseFrequencyOperatorTest.class,
	SingleBaseOperatorTest.class, 
	SingleBaseUniformOperatorTest.class,


})
public class AAllTestsHaplotyesOperator {

}
