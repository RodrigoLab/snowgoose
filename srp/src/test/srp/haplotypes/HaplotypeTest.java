package test.srp.haplotypes;

import static org.junit.Assert.*;

import java.lang.reflect.Method;
import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import srp.evolution.haplotypes.old.OldHaplotype;
import srp.haplotypes.Haplotype;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

public class HaplotypeTest {
	
	@Rule public ExpectedException thrown= ExpectedException.none();
	
	public static final Class<Haplotype> haplotypeClass = Haplotype.class;
	public static final String className = haplotypeClass.getName();
		
	private Haplotype haplotype;

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		haplotype = new Haplotype("ACGT");
		
	}

	@After
	public void tearDown() throws Exception {
	}
	
	
	@Test
	public void testOverride() throws Exception {
		
		//Final method cannot be override
		assertFalse( isMethodOverriden("setState", int.class, int.class));

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
	
	@Test
	public void testHaplotype() throws Exception {
		assertEquals(4, haplotype.getLength());
		assertEquals("ACGT", haplotype.getSequenceString());
		
		haplotype = new Haplotype(new Sequence("AACCGGTT"));
		assertEquals("AACCGGTT", haplotype.getSequenceString());
		assertEquals(8, haplotype.getLength());
		
		haplotype = new Haplotype(new Taxon("testTaxonID"), "TTGGCCAACCGGTT");
		assertEquals("TTGGCCAACCGGTT", haplotype.getSequenceString());
		assertEquals(14, haplotype.getLength());
		assertEquals("testTaxonID", haplotype.getTaxon().getId());
	}

	@Test
	public void testSetGetChar() throws Exception {
		assertEquals('A', haplotype.getChar(0));
		assertEquals('C', haplotype.getChar(1));
		assertEquals('G', haplotype.getChar(2));
		assertEquals('T', haplotype.getChar(3));

		haplotype.setCharAt(0, 'T');
		haplotype.setCharAt(1, (int) 'G');
		haplotype.setCharAt(2, 'C');
		haplotype.setCharAt(3, (int) 'A');
		
		assertEquals('A', haplotype.getChar(3));
		assertEquals('C', haplotype.getChar(2));
		assertEquals('G', haplotype.getChar(1));
		assertEquals('T', haplotype.getChar(0));
		
		haplotype = new Haplotype("ACGTACGT");
		char[] charArrays = new char[10];
		Arrays.fill(charArrays, 'X');
		haplotype.getChars(1, 6, charArrays, 2);
		char[] expecteds = new char[]{'X', 'X', 'C', 'G', 'T', 'A', 'C', 'X', 'X', 'X'};
		assertArrayEquals(expecteds, charArrays);
		

	}

	public static boolean isMethodOverriden(String methodName, Class<?>... parameterTypes ) throws NoSuchMethodException, SecurityException{
		
	//		Class<?> declaringClass = haplotypeClass.getMethod(methodName, parameterTypes ).getDeclaringClass();
	//		return declaringClass.getName().equals(className);
			Class<?> declaringClass = OldHaplotype.class.getMethod(methodName, parameterTypes ).getDeclaringClass();
			return declaringClass.getName().equals(OldHaplotype.class.getName());
		}

	public static boolean isMethodOverrridenGeneric(Method myMethod) {
		Class<?> declaringClass = myMethod.getDeclaringClass();
		System.out.println(declaringClass.getCanonicalName());
		if (declaringClass.equals(Object.class)) {
			return false;
		}
		try {
			declaringClass.getSuperclass().getMethod(myMethod.getName(),
					myMethod.getParameterTypes());
			return true;
		} catch (NoSuchMethodException e) {
			return false;
		}
	}

	@Test
	public void testGetState() throws Exception {
		for (int i = 0; i < 4; i++) {
			assertEquals(i, haplotype.getState(i));
		}
		
	}

	@Test
	public void testGuessDataType() throws Exception {
		DataType guessDataType = haplotype.guessDataType();
		assertEquals(Nucleotides.INSTANCE.getName(), guessDataType.getName());
	}

	@Test
	public void testSetSequenceString() throws Exception {
		haplotype.setSequenceString("TTGG");
		assertEquals("TTGG", haplotype.getSequenceString());;
		thrown.expect(IllegalArgumentException.class);
		haplotype.setSequenceString("TTGGCC");
		
		
	}

}
