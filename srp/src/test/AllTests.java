package test;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;


@RunWith(Suite.class)
@SuiteClasses({ 
//	test.dr.ext.TreeLikelihoodExtTest.class, //slow
	test.srp.haplotypes.AAllTestsHaplotypes.class,
	test.srp.operator.haplotypes.AAllTestsHaplotyesOperator.class,
	test.srp.likelihood.haplotypes.AllTestsLikelihood.class,
	test.srp.shortreads.AAllTestsShortreads.class,
	test.srp.spectrum.AAllTestsSpectrum.class,
	test.srp.likelihood.spectrum.AllTestsSpectrumLikelihood.class,
	test.srp.operator.spectrum.AllTestsSpectrumOperator.class,

	})
public class AllTests {

}
