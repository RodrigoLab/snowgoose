package test.srp.haplotypes.operator;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ 
	HaplotypeRecombinationOperatorTest.class,
	HaplotypeSwapSectionOperatorTest.class,
	SingleBaseOperatorTest.class, 
	SwapBasesEmpiricalOperatorTest.class,
	SwapBasesMultiOperatorTest.class,
	SwapBaseUniformOperatorTest.class,

})
public class AllTestsHaplotyesOperator {

}
