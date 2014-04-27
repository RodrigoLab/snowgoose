package test.srp.dbgraph;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ 

	CompatibleNodeTest.class,
	CompatibleSetsTest.class,
	DeBruijnGraphLikelihoodTest.class,
	DeBruijnGraphTest.class,
	DeBruijnImporterTest.class
	
	})
public class AllTestDBGraph {

}
