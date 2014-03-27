package test.srp.spectrum;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ 
	AbstractSpectraTest.class,
	AbstractSpectrumAlignmentModelTest.class,
	CategorySpectraParameterTest.class,
	CategorySpectrumAlignmentModelTest.class,
	CategorySpectrumTest.class,
	SpectraParameterTest.class,
	SpectrumAlignmentModelTest.class,
	SpectrumAlignmentUtilsTest.class,
	SpectrumTest.class,
})
public class AAllTestsSpectrum {

}
