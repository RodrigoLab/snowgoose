package test.srp.haplotypes.old;

import static org.junit.Assert.assertTrue;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.evolution.haplotypes.old.OldHaplotype;

public class OldHaplotypeTest {
	
	public static final Class<OldHaplotype> oldHaplotypeClass = OldHaplotype.class;
	public static final String className = oldHaplotypeClass.getName();

	private OldHaplotype haplotype;
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		haplotype = new OldHaplotype("ACTG");
		
	}

	@After
	public void tearDown() throws Exception {
	}
	
	
	@Test
	public void testOverride() throws Exception {
		
		//Final method cannot be override
		assertTrue(! isMethodOverriden("setState", int.class, int.class));
		
		assertTrue(isMethodOverriden("getLength"));
		assertTrue(isMethodOverriden("getSequenceString"));
		assertTrue(isMethodOverriden("getChar", int.class));
		assertTrue(isMethodOverriden("getState", int.class));
		assertTrue(isMethodOverriden("getChars", int.class, int.class, char[].class, int.class));
		assertTrue(isMethodOverriden("guessDataType"));
		assertTrue(isMethodOverriden("setSequenceString", String.class));
		assertTrue(isMethodOverriden("appendSequenceString", String.class));
		assertTrue(isMethodOverriden("insertSequenceString", int.class, String.class));

//		Method[] methods = haplotypeClass.getMethods();
//		for (int i = 0; i < methods.length; i++) {
//			System.out.println(methods[i].getName());
//			System.out.println(methods[i].getDeclaringClass().getName());
//			System.out.println();
//		}

		
	}
	
	public static boolean isMethodOverriden(String methodName, Class<?>... parameterTypes ) throws NoSuchMethodException, SecurityException{
	
		Class<?> declaringClass = oldHaplotypeClass.getMethod(methodName, parameterTypes ).getDeclaringClass();
		return declaringClass.getName().equals(className);
//		Class<?> declaringClass = OldHaplotype.class.getMethod(methodName, parameterTypes ).getDeclaringClass();
//		return declaringClass.getName().equals(OldHaplotype.class.getName());
	}

}
