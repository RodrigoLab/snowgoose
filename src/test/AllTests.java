package test;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

import test.alignment.AAllTestsAlignment;


@RunWith(Suite.class)
@SuiteClasses({ 
	test.likelihood.AllTestsLikelihood.class,
	test.alignment.AAllTestsAlignment.class
	})
public class AllTests {

}
