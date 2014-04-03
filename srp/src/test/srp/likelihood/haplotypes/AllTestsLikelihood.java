package test.srp.likelihood.haplotypes;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

import test.srp.likelihood.LikelihoodScalerTest;

@RunWith(Suite.class)
@SuiteClasses({ 

	LikelihoodScalerTest.class,
	LikelihoodUtilsTest.class, 
	ShortReadLikelihoodTest.class,
	ShortReadLikelihoodWithOperatorTest.class
})
public class AllTestsLikelihood {

}
