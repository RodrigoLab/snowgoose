package test;

import static org.junit.Assert.assertTrue;
import dr.inference.markovchain.MarkovChain;

public class TestUtils {

	public static final double UNITTEST_THRESHOLD = 1e-8;
	public static final double THRESHOLD = MarkovChain.EVALUATION_TEST_THRESHOLD;
	public static void assertExpectationRange(double mean, double value, double error) {
        double upper = value + error;
        double lower = value - error;

        assertTrue("Expected value is: " + value + " but got " + mean + " +/- " + error,
                upper > mean && lower < mean);
    }
	
	public static void assertExpectationStderr(String name, double mean, double value, double stderr) {
        double upper = value + 2*stderr;
        double lower = value - 2*stderr;

        assertTrue("Expected " + name + " is " + value + " but got " + mean + " +/- 2*" + stderr,
                upper > mean && lower < mean);
    }
}
