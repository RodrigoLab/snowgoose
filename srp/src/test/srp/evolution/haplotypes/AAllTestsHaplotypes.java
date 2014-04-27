package test.srp.evolution.haplotypes;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

import srp.evolution.haplotypes.old.OldHaplotype;
import test.srp.evolution.haplotypes.old.OldHaplotypeModelTest;

@RunWith(Suite.class)
@SuiteClasses({ 
	AbstractHaplotypeModelTest.class,
 
	AlignmentUtilsTest.class,

	SPSDistTest.class, 
//	ShortReadTest.class
	
	OldHaplotypeModelTest.class,
	OldHaplotype.class,
	})
public class AAllTestsHaplotypes {

}
