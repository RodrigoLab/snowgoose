package test.srp.haplotypes.operator;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ 
	AbstractSwapBasesOperatorTest.class,
	
	HaplotypeRecombinationOperatorTest.class,
	HaplotypeSwapSectionOperatorTest.class,
	
	SingleBaseFrequencyOperatorTest.class,
	SingleBaseOperatorTest.class, 
	SingleBaseUniformOperatorTest.class,

	SwapBasesEmpiricalOperatorTest.class,
	SwapBasesMultiOperatorTest.class,
	SwapBasesUniformOperatorTest.class,

})
public class AllTestsHaplotyesOperator {

}
