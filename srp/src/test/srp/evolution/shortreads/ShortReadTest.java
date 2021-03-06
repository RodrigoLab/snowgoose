package test.srp.evolution.shortreads;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.evolution.shortreads.ShortRead;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

public class ShortReadTest {

	private ShortRead srp;
	private static Taxon taxon;

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		
		taxon = new Taxon("id");
		String seqs = "........ACGTACGT****ACGTACGT....";
		Sequence seq = new Sequence(taxon, seqs);
		srp = new ShortRead(seq);
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testGetFullSrpCharAt() {
		assertEquals('.', srp.getFullSrpCharAt(0));
		assertEquals('.', srp.getFullSrpCharAt(7));
		assertEquals('A', srp.getFullSrpCharAt(8));
		assertEquals('C', srp.getFullSrpCharAt(9));
		assertEquals('G', srp.getFullSrpCharAt(10));
		assertEquals('T', srp.getFullSrpCharAt(11));
	}

	@Test
	public void testGetSrpString() {
		String expected = "........ACGTACGT----ACGTACGT....";
		assertEquals(expected, srp.getFullSrp());
		expected = "ACGTACGT----ACGTACGT";
		assertEquals(expected, srp.getFragmentSrp());
		
		expected = "id";
		assertEquals(expected, srp.getName());
		
	}
	
	@Test
	public void testPosition(){
//		s = "........ACGTACGT****ACGTACGT....";
		assertEquals(8, srp.getStart());
		assertEquals(28, srp.getEnd());
	}
	
	@Test
	public void testIsValid(){
		
		String s = "...***...";
		ShortRead srTemp = new ShortRead(new Sequence(taxon, s));
		assertTrue(srTemp.getIsValid());
		assertEquals("---", srTemp.getFragmentSrp());
		
		s = "AA";
		srTemp = new ShortRead(new Sequence(taxon, s));
		assertTrue(srTemp.getIsValid());
		assertEquals("AA", srTemp.getFragmentSrp());
		
		s = "...AAA";
		srTemp = new ShortRead(new Sequence(taxon, s));
		assertTrue(srTemp.getIsValid());
		assertEquals("AAA", srTemp.getFragmentSrp());
		
		s = "AAAA...";
		srTemp = new ShortRead(new Sequence(taxon, s));
		assertTrue(srTemp.getIsValid());
		assertEquals("AAAA", srTemp.getFragmentSrp());
		
		s = "..ZZ..";
		srTemp = new ShortRead(new Sequence(taxon, s));
		assertFalse(srTemp.getIsValid());
		assertEquals("..ZZ..", srTemp.getFullSrp());
		
	}

}
