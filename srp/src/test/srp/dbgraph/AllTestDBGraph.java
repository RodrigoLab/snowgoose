package test.srp.dbgraph;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

import srp.evolution.haplotypes.old.OldHaplotype;
import test.srp.haplotypes.old.OldHaplotypeModelTest;
import test.srp.shortreads.AlignmentMappingTest;
import test.srp.shortreads.ShortReadTest;

@RunWith(Suite.class)
@SuiteClasses({ 

	CompatibleNodeTest.class,
	CompatibleSetsTest.class,
	DeBruijnGraphTest.class,
	DeBruijnImporterTest.class
	
	})
public class AllTestDBGraph {

}
