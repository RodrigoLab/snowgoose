package test.srp.spectrum.operator;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ 
	ColumnSpectrumDeltaExchangeOperatorTest.class,
	MultiSpectrumDeltaExchangeOperatorTest.class,
	RecombinationSpectrumOperatorTest.class,
	SingleSpectrumDeltaExchangeOperatorTest.class,
})
public class AllTestsSpectrumOperator {

}
