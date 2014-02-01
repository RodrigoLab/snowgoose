package test.srp.haplotypes.likelihood;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ 

	LikelihoodScalerTest.class,
	LikelihoodUtilsTest.class, 
	ShortReadLikelihoodTest.class,
	ShortReadLikelihoodWithOperatorTest.class
})
public class AllTestsLikelihood {

}
