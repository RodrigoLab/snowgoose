package test.alignment;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({ AlignmentMappingTest.class, AlignmentUtilsTest.class,
		HaplotypesTest.class, ShortReadTest.class })
public class AAllTestsAlignment {

}
