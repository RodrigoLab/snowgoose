package test;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;


@RunWith(Suite.class)
@SuiteClasses({ 
//	test.dr.ext.TreeLikelihoodExtTest.class, //slow
	test.srp.haplotypes.AAllTestsHaplotypes.class,
	test.srp.haplotypes.operator.AAllTestsHaplotyesOperator.class,
	test.srp.likelihood.AllTestsLikelihood.class,
	
	})
public class AllTests {

}
