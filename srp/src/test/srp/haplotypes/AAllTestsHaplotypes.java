package test.srp.haplotypes;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ 
	AbstractHaplotypeModelTest.class,
	AlignmentMappingTest.class, 
	AlignmentUtilsTest.class,
	HaplotypeModelTest.class,
	HaplotypeModelUtilsTest.class, 
	ShortReadTest.class })
public class AAllTestsHaplotypes {

}
