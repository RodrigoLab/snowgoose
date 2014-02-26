package test;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;


@RunWith(Suite.class)
@SuiteClasses({ 
//	test.dr.ext.TreeLikelihoodExtTest.class, //slow
	test.srp.haplotypes.AAllTestsHaplotypes.class,
	test.srp.haplotypes.operator.AAllTestsHaplotyesOperator.class,
	test.srp.haplotypes.likelihood.AllTestsLikelihood.class,
	test.srp.shortreads.AAllTestsShortreads.class,
	test.srp.spectrum.AAllTestsSpectrum.class,
	test.srp.spectrum.likelihood.AllTestsSpectrumLikelihood.class,
	test.srp.spectrum.operator.AllTestsSpectrumOperator.class,

	})
public class AllTests {

}
