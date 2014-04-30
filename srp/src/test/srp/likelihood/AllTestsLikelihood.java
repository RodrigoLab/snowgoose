package test.srp.likelihood;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

import test.srp.likelihood.haplotypes.AllTestsLikelihoodHaplotypes;
import test.srp.likelihood.haplotypes.LikelihoodUtilsTest;
import test.srp.likelihood.haplotypes.ShortReadsHaplotypeLikelihoodTest;
import test.srp.likelihood.spectrum.AllTestsSpectrumLikelihood;

@RunWith(Suite.class)
@SuiteClasses({ 

	LikelihoodScalerTest.class,
	AllTestsLikelihoodHaplotypes.class,
	AllTestsSpectrumLikelihood.class,
})
public class AllTestsLikelihood {

}
