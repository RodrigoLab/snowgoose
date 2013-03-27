package test;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;


@RunWith(Suite.class)
@SuiteClasses({ 
	test.likelihood.AllTestsLikelihood.class,
	test.haplotypes.AAllTestsAlignment.class
	})
public class AllTests {

}
