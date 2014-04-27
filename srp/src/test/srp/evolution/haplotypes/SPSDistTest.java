package test.srp.evolution.haplotypes;


import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.evolution.haplotypes.SPSDist;
import srp.evolution.haplotypes.old.OldHaplotypeModel;
import srp.evolution.haplotypes.old.OldHaplotypeModelUtils;


public class SPSDistTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}



	@Test
	public void testCalculateSpsDoubleSame() {
	
		String[] seq1 = new String[]{
				"AACCTTGG",
				"ACCTTGGA",
				"CCTTGGAA"
		};
		String[] seq2 = new String[]{
				"AACCTTGG",
				"ACCTTGGA",
				"CCTTGGAA"
		};
		OldHaplotypeModel h1 = OldHaplotypeModelUtils.createHaplotypeModel(seq1);
		OldHaplotypeModel h2 = OldHaplotypeModelUtils.createHaplotypeModel(seq2);
		int sps = SPSDist.calculeteSPS(h1, h2);
		int expected = 4+8 + 4+4 + 4+8;
		assertEquals(expected, sps, 0);
		
		int[][] spsArray = SPSDist.calculeteSPSArray(h1, h2);
		int[][] expectedArray = new int[][]{
					{ 0, 4, 8 }, 
					{ 4, 0, 4 },
					{ 8, 4, 0 }		
				};
		
		seq2 = new String[]{
					"AACTCTGA", //3+3+6
					"AAGTCCTA", //6+5+6
					"AACAATGA"	//3+4+7
				};
		h2 = OldHaplotypeModelUtils.createHaplotypeModel(seq2);
		sps = SPSDist.calculeteSPS(h1, h2);
		expected = 43;
		assertEquals(expected, sps);

		spsArray = SPSDist.calculeteSPSArray(h1, h2);
		expectedArray = new int[][] { 
						{ 3, 6, 3 }, 
						{ 3, 5, 4 }, 
						{ 6, 6, 7 } 
				};
		assertArrayEquals(expectedArray, spsArray);
	}

}
