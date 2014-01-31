package test.srp.spectrum.likelihood;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

import test.srp.spectrum.treelikelihood.SpectrumTreeLikelihoodTest;

@RunWith(Suite.class)
@SuiteClasses({ 
	ShortReadsSpectrumLikelihoodTest.class,
	ShortReadsSpectrumLikelihoodTimeTrialTest.class,
	SpectrumTreeLikelihoodTest.class
	})
public class AllTestsSpectrumLikelihood {

}
