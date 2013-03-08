package test.srp;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;


@RunWith(Suite.class)
@SuiteClasses({ 
	test.srp.likelihood.AllTestsLikelihood.class,
	test.srp.haplotypes.AAllTestsAlignment.class
	})
public class AllTests {

}
