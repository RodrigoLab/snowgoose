package test.srp.likelihood;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ 

	LikelihoodScalerTest.class,
	LikelihoodUtilsTest.class, 
	ShortReadLikelihoodTest.class
})
public class AllTestsLikelihood {

}
