package test.srp.haplotypes.operator;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ 
	HaplotypeRecombinationOperatorTest.class,
	HaplotypeSwapSectionOperatorTest.class,
	SwapBaseOperatorTest.class, 
	SwapMultiBasesOperatorTest.class,
	UniformSwapBaseOperatorTest.class,

})
public class AllTestsHaplotyesOperator {

}
