package test.srp.likelihood.spectrum;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

import test.srp.likelihood.spectrum.treelikelihood.SpectrumTreeLikelihoodTest;
import test.srp.likelihood.stateLikelihood.StateLikelihoodTest;

@RunWith(Suite.class)
@SuiteClasses({ 
	ShortReadsSpectrumLikelihoodTest.class,
//	ShortReadsSpectrumLikelihoodTimeTrialTest.class,

	StateLikelihoodTest.class,

	SpectrumTreeLikelihoodTest.class
	})
public class AllTestsSpectrumLikelihood {

}
