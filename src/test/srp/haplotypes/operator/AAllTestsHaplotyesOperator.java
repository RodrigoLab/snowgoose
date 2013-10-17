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
	
	SingleBaseFrequencyOperatorTest.class,
	SingleBaseOperatorTest.class, 
	SingleBaseUniformOperatorTest.class,

	MultiBasesEmpiricalOperatorTest.class,
	MultiBasesOperatorTest.class,
	MultiBasesUniformOperatorTest.class,

})
public class AAllTestsHaplotyesOperator {

}
