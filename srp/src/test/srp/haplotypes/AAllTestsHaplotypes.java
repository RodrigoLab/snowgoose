package test.srp.haplotypes;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

import test.srp.haplotypes.old.OldHaplotypeModelTest;
import test.srp.shortreads.AlignmentMappingTest;
import test.srp.shortreads.ShortReadTest;

@RunWith(Suite.class)
@SuiteClasses({ 
	AbstractHaplotypeModelTest.class,
 
	AlignmentUtilsTest.class,
//	OldHaplotypeModelTest.class,
	SPSDistTest.class, 
//	ShortReadTest.class 
	})
public class AAllTestsHaplotypes {

}
