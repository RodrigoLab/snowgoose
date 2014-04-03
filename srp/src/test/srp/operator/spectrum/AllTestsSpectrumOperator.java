package test.srp.operator.spectrum;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ 
	DeltaExchangeColumnSpectrumOperatorTest.class,
	DeltaExchangeMultiSpectrumOperatorTest.class,
	DeltaExchangeSingleSpectrumOperatorTest.class,
	DirichletAlphaSpectrumOperatorTest.class,
	DirichletSpectrumOperatorTest.class,
	RecombineSectionSpectrumOperatorTest.class,
	SwapSingleSpectrumOperatorTest.class
	
})
public class AllTestsSpectrumOperator {

}
