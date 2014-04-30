package test;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;


@RunWith(Suite.class)
@SuiteClasses({ 
//	test.dr.ext.TreeLikelihoodExtTest.class, //slow
	test.srp.evolution.haplotypes.AAllTestsHaplotypes.class,
	test.srp.evolution.shortreads.AAllTestsShortreads.class,
	test.srp.evolution.spectrum.AAllTestsSpectrum.class,
	
	test.srp.likelihood.AllTestsLikelihood.class,
//	test.srp.likelihood.spectrum.AllTestsSpectrumLikelihood.class,

	test.srp.operator.haplotypes.old.AAllTestsHaplotyesOperator.class,
	test.srp.operator.spectrum.AllTestsSpectrumOperator.class,

	})
public class AllTests {

}
