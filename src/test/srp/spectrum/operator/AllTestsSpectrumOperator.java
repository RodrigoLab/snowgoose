package test.srp.spectrum.operator;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ 
	DeltaExchangeColumnSpectrumOperatorTest.class,
	DeltaExchangeMultiSpectrumOperatorTest.class,
	RecombineSectionSpectrumOperatorTest.class,
	DeltaExchangeSingleSpectrumOperatorTest.class,
})
public class AllTestsSpectrumOperator {

}
