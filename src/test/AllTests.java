package test;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;


@RunWith(Suite.class)
@SuiteClasses({ 
//	test.dr.ext.TreeLikelihoodExtTest.class,
	test.haplotypes.AAllTestsAlignment.class,
//	test.haplotypes.operator.AlignmentSwapBaseOperatorTest.class,
	test.likelihood.AllTestsLikelihood.class,
	
	})
public class AllTests {

}
