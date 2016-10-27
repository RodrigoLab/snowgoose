package test.srp.evolution.haplotypes;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

import srp.evolution.haplotypes.old.OldHaplotype;
import test.srp.evolution.haplotypes.old.OldHaplotypeModelTest;
import test.srp.evolution.haplotypes.old.OldHaplotypeTest;

@RunWith(Suite.class)
@SuiteClasses({ 
	AbstractHaplotypeModelTest.class,
 
	AlignmentUtilsTest.class,

	HaplotypeModelTest.class,
	HaplotypeTest.class,
	SPSDistTest.class, 
//	ShortReadTest.class
	
	OldHaplotypeModelTest.class,
	OldHaplotypeTest.class,
	})
public class AAllTestsHaplotypes {

}
