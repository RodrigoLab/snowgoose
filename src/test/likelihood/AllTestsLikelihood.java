package test.likelihood;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ 
	LikelihoodCalculationTest.class, 
	LikelihoodScalerTest.class,
	LikelihoodUtilsTest.class, 
	ShortReadLikelihoodTest.class
})
public class AllTestsLikelihood {

}
