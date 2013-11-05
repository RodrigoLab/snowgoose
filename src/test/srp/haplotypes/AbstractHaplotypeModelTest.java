package test.srp.haplotypes;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.hamcrest.CoreMatchers;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import srp.core.DataImporter;
import srp.dr.evolution.alignment.ShortReadFragmentsAlignment;
import srp.dr.evolution.datatype.ShortReads;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.Haplotype;
import srp.haplotypes.HaplotypeModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxon;


public class AbstractHaplotypeModelTest {

	private static AlignmentMapping aMap;
	private static Alignment srpAlignment;
	
	private static Taxon[] expectedTaxons;
	private static String[] expectedSequences;
	private static ArrayList<Haplotype> expectedList;
	private static HashMap<Character, Integer> charToState;
	
	private HaplotypeModel haplotypeModel;
	private HaplotypeModel haplotypeModelRandom;

	
	
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {

		String dir = System.getProperty("user.dir")+File.separatorChar+"unittest"+File.separator;
		srpAlignment = DataImporter.importShortReads(dir, "HaplotypeModelTest_10_srp.fasta");
		aMap = new AlignmentMapping(srpAlignment);
		

		expectedSequences = new String[]{ 
				".............................................AAGCTCCAACACA*GAGAGCAACTCGTAGGGCGTTGTTCAA..............",
				"..............CCGCGCC*TGTCAGGGCTCAAAACCC*GA*AAAGC...................................................",
				"..............................CTCAAAATCCTGAGAAAGCTCCAACGCACGAGAGCAACT...............................",
				".................CGCC*TGTCAGGGCTCAAAATCCTGAGAAAGCTCC................................................",
				".........................TAGGGCTCAAAATCCTGAGAAAGCTCCAACACACGAG......................................",
				"CCAGG...............................................................................................",
				"...............................................GCTCCAACACACGAGAGC*ACTCGTAGGGCGTTGTTCAATTC...........",
				".......GTCCAGTCCGCGCCTTGTCAGGGCTCG..................................................................",
				"........................................................CACGAGAGCAACACGTA*GGCGTTG*TCAATTCTAGTTCT....",
				"..........CAGTCCGCGCCTTGTCA*GGCTCAAAAT*CTG..........................................................",
		};
		expectedTaxons = new Taxon[10];
//		char[][] expectedMatrix = new char[matrixS.length][matrixS[0].length()];
		expectedList = new ArrayList<Haplotype>();
		for (int i = 0; i < expectedSequences.length; i++) {
			expectedTaxons[i] = new Taxon("r"+(i+1)+".1");
			Haplotype h = new Haplotype(expectedTaxons[i], expectedSequences[i]);
			expectedList.add(h);
		}
		
		charToState = new HashMap<Character, Integer>();
		charToState.put('A', 0);
		charToState.put('C', 1);
		charToState.put('G', 2);
		charToState.put('T', 3);
		charToState.put('-', 17);
		charToState.put('*', 17);
		charToState.put('.', 17);
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		haplotypeModel = new HaplotypeModel(aMap, srpAlignment);
		haplotypeModelRandom = new HaplotypeModel(aMap, 5);
	}

	@After
	public void tearDown() throws Exception {
	}


	@Test
	public void testHaplotypes() throws Exception {
		
		for (int i = 0; i < expectedList.size(); i++) {
		
			assertEquals(expectedSequences[i], haplotypeModel.getAlignedSequenceString(i));
			assertEquals(expectedSequences[i], haplotypeModel.getUnalignedSequenceString(i));
			assertEquals(expectedSequences[i], haplotypeModel.getHaplotype(i).getSequenceString());
			assertEquals(expectedTaxons[i], haplotypeModel.getHaplotype(i).getTaxon());
			
			assertEquals(haplotypeModel.getSequence(i), haplotypeModel.getHaplotype(i));
//			expectedMatrix[i] = matrixS[i].toCharArray();
			for (int j = 0; j < expectedSequences[i].length(); j++) {
				char expected = expectedSequences[i].charAt(j); 
				assertEquals(expected,haplotypeModel.getHaplotypeCharAt(i,j));
			}
		}
//		assertArrayEquals(expectedMatrix, haplotypeModel.getCharMatrix());
		
	}
	
	
	@Test
	public void testSequenceList() throws Exception {
		// getSequence(int)
		// in testHaplotypes()
		
		// getSequenceCount()
		// in testInit()
		
		// getSequenceAttribute(int, String)
		// setSequenceAttribute(int, String, Object)
		for (int i = 0; i < haplotypeModel.getSequenceCount(); i++) {
			assertNull(haplotypeModel.getSequenceAttribute(i, "NULL"));
			haplotypeModel.setSequenceAttribute(i, "NAME", expectedTaxons[i].getId());
		}
		for (int i = 0; i < haplotypeModel.getSequenceCount(); i++) {
			assertEquals(("r"+(i+1)+".1"), haplotypeModel.getSequenceAttribute(i, "NAME"));
			assertNull(haplotypeModel.getSequenceAttribute(i, "NULL"));
		}
		
	}
	@Test
	public void testTaxonList() throws Exception {

		// getTaxonCount()
		assertEquals(10, haplotypeModel.getTaxonCount());
		assertEquals(5, haplotypeModelRandom.getTaxonCount());

		// getTaxon(int)
		// getTaxonId(int)
		// getTaxonIndex(String)
		// getTaxonIndex(Taxon)
		
		for (int i = 0; i < expectedTaxons.length; i++) {
			String expectedID = expectedTaxons[i].getId();
			assertEquals(expectedTaxons[i], haplotypeModel.getTaxon(i));
			assertEquals(expectedID, haplotypeModel.getTaxonId(i));
			assertEquals(i, haplotypeModel.getTaxonIndex(expectedID));
			assertEquals(i, haplotypeModel.getTaxonIndex(haplotypeModel.getTaxon(i)));
			
		}
		assertEquals(-1, haplotypeModel.getTaxonIndex("EMPTY"));
		assertEquals(-1, haplotypeModel.getTaxonIndex(new Taxon("E")));
		assertEquals(-1, haplotypeModel.getTaxonIndex(expectedTaxons[0]));
		
		List<Taxon> expectedList = new ArrayList<Taxon>();
		for (int i = 0; i < haplotypeModelRandom.getTaxonCount(); i++) {
			String expectedID = HaplotypeModel.TAXON_PREFIX+i;
			Taxon expectedTaxon = new Taxon(expectedID);
			expectedList.add(expectedTaxon);
			assertEquals(expectedTaxon, haplotypeModelRandom.getTaxon(i));
			assertEquals(expectedID, haplotypeModelRandom.getTaxonId(i));
			assertEquals(i, haplotypeModelRandom.getTaxonIndex(expectedID));
			assertEquals(i, haplotypeModelRandom.getTaxonIndex(haplotypeModelRandom.getTaxon(i)));
		}
		assertEquals(-1, haplotypeModelRandom.getTaxonIndex("EMPTY"));
		assertEquals(-1, haplotypeModelRandom.getTaxonIndex(expectedTaxons[0]));
		
		
		// asList()		
		assertEquals(Arrays.asList(expectedTaxons), haplotypeModel.asList());
		assertEquals(expectedList, haplotypeModelRandom.asList());

		
		// getTaxonAttribute(int, String)
		for (int i = 0; i < haplotypeModelRandom.getTaxonCount(); i++) {
			Taxon taxon = haplotypeModelRandom.getTaxon(i);
			taxon.setAttribute("NAME", expectedTaxons[i].getId());
			taxon.setAttribute("INDEX", i);
		}
		for (int i = 0; i < haplotypeModelRandom.getTaxonCount(); i++) {
			assertEquals(("r"+(i+1)+".1"), haplotypeModelRandom.getTaxonAttribute(i, "NAME"));
			assertEquals(i, haplotypeModelRandom.getTaxonAttribute(i, "INDEX"));
			assertNull(haplotypeModelRandom.getTaxonAttribute(i, "NULL"));
			
		}
		
		SimpleAlignment alignmentNoTaxon = new SimpleAlignment();
		alignmentNoTaxon.setDataType(ShortReads.INSTANCE);
		Haplotype h;
		for (int i = 0; i < 5; i++) {
			if (i == 0) {
				h = AbstractHaplotypeModelTest.expectedList.get(0);
			} else {
				h = new Haplotype(expectedSequences[i]);
			}

			alignmentNoTaxon.addSequence(h);

		}
		haplotypeModelRandom = new HaplotypeModel(aMap, alignmentNoTaxon);
		assertEquals(5, haplotypeModelRandom.getTaxonCount());
		for (int i = 1; i < haplotypeModelRandom.getTaxonCount(); i++) {
			
			assertNull(haplotypeModelRandom.getTaxon(i));
			
			assertNull(haplotypeModelRandom.getTaxonAttribute(i, "INDEX"));
			haplotypeModelRandom.setSequenceAttribute(i, "INDEX", i+i);
			
			assertEquals(i+i, haplotypeModelRandom.getSequenceAttribute(i,"INDEX"));
			assertEquals(i+i, haplotypeModelRandom.getTaxonAttribute(i, "INDEX"));
			
		}
		
		for (int i = 0; i < haplotypeModelRandom.getTaxonCount(); i++) {

			try {
				haplotypeModelRandom.getTaxonId(i);
				
			} catch (IllegalArgumentException e) {
		        assertThat(e.getMessage(), CoreMatchers.startsWith("Illegal taxon index"));
		        assertThat(e.getMessage(), CoreMatchers.containsString(""+i));
			}
		
		}
	}
	
	@Rule
 	public ExpectedException exception = ExpectedException.none();

//	
	@Test
	public void testSiteList() throws Exception {
		
		// getPatternIndex(int)
		Random rand = new Random();
		for (int i = 0; i < 10; i++) {
			int expected = rand.nextInt(100);
			assertEquals(expected, haplotypeModel.getPatternIndex(expected));
			assertEquals(expected, haplotypeModelRandom.getPatternIndex(expected));
		}
		
		// getSiteCount()
		assertEquals(100, haplotypeModel.getSiteCount());
		assertEquals(100, haplotypeModelRandom.getSiteCount());
		
		// getSitePattern(int)
		int[] expectedArray = new int[10];
		for (int site = 0; site < expectedSequences[0].length(); site++) {
			for (int i = 0; i < expectedArray.length; i++) {
				expectedArray[i] = charToState.get(expectedSequences[i].charAt(site));
			}
			assertArrayEquals(expectedArray, haplotypeModel.getSitePattern(site));
		}
		Arrays.fill(expectedArray, 17);
		assertArrayEquals(expectedArray, haplotypeModel.getSitePattern(5));
		assertArrayEquals(expectedArray, haplotypeModel.getSitePattern(101));
		
		

		
		// getState(int, int)
		for (int h = 0; h < haplotypeModelRandom.getHaplotypeCount(); h++) {
			String hap = haplotypeModelRandom.getHaplotypeString(h);
			assertEquals(hap.length(), haplotypeModelRandom.getHaplotypeLength());
			assertEquals(-1, hap.indexOf('.'));
			for (int i = 0; i < hap.length(); i++) {
				int expected = charToState.get(hap.charAt(i));
				assertEquals(expected, haplotypeModelRandom.getState(h, i));
				assertEquals(expected, haplotypeModelRandom.getPatternState(h, i));
			}
			int expected = 17; 
			assertEquals(expected, haplotypeModelRandom.getState(h, hap.length()));
			assertEquals(expected, haplotypeModelRandom.getPatternState(h, hap.length()));
		}
		
		for (int h = 0; h < haplotypeModel.getHaplotypeCount(); h++) {
			String hap = expectedSequences[h];
			assertEquals(hap.length(), haplotypeModel.getHaplotypeLength());
			for (int i = 0; i < hap.length(); i++) {
				int expected = charToState.get(hap.charAt(i));
				assertEquals(expected, haplotypeModel.getState(h, i));
				assertEquals(expected, haplotypeModel.getPatternState(h, i));
			}
			int expected = 17;
			assertEquals(expected, haplotypeModel.getState(h, hap.length()));
			assertEquals(expected, haplotypeModel.getPatternState(h, hap.length()));
		}

	}

	@Test
	public void testPatternList() throws Exception {
		// getDataType()
		assertEquals(Nucleotides.INSTANCE, haplotypeModel.getDataType());
		assertEquals(Nucleotides.class, haplotypeModel.getDataType().getClass());
		assertEquals(Nucleotides.DESCRIPTION, haplotypeModel.getDataType().getDescription());

		// getPattern(int)
		// in testSiteList()
		
		// getPatternCount()
		assertEquals(100, haplotypeModel.getPatternCount());
		assertEquals(100, haplotypeModelRandom.getPatternCount());
		
		// getPatternLength()
		assertEquals(10, haplotypeModel.getPatternLength());
		assertEquals(5, haplotypeModelRandom.getPatternLength());
		
		// getPatternState(int, int)
		// in testSiteList()
		
		
		// getPatternWeight(int)
		for (int i = 0; i < 100; i++) {
			assertEquals(1, haplotypeModel.getPatternWeight(i), 0);
			assertEquals(1, haplotypeModelRandom.getPatternWeight(i), 0);
		}
		
		// getPatternWeights()
		double[] expectedArray = new double[haplotypeModel.getPatternCount()];
		Arrays.fill(expectedArray, 1.0);
		assertArrayEquals(expectedArray, haplotypeModel.getPatternWeights(), 0);
		assertArrayEquals(expectedArray, haplotypeModelRandom.getPatternWeights(), 0);
		
		// getStateCount()
		assertEquals(4, haplotypeModel.getStateCount());
		assertEquals(4, haplotypeModelRandom.getStateCount());
		
		// getStateFrequencies()
		double[] expectedFreq = {0.2755418, 0.3034056, 0.2414861, 0.1795666};
		assertArrayEquals(expectedFreq, haplotypeModel.getStateFrequencies(), 1e-7);
		
		double[] stateFrequencies = haplotypeModelRandom.getStateFrequencies();
		StringBuilder sb = new StringBuilder();
		
		for (int i = 0; i < haplotypeModelRandom.getHaplotypeCount(); i++) {
			sb.append(haplotypeModelRandom.getHaplotypeString(i));
		}
		String temp = sb.toString();
		double[] count = new double[4];
		Arrays.fill(count, 0);
		for (int i = 0; i < temp.length(); i++) {
			char c = temp.charAt(i);
			switch (c) {
			case 'A':
				count[0]++;
				break;
			case 'C':
				count[1]++;
				break;
			case 'G':
				count[2]++;
				break;
			case 'T':
				count[3]++;
				break;
			default:
				break;
			}
		}
		double total = StatUtils.sum(count);
		for (int i = 0; i < count.length; i++) {
			count[i] /= total;
		}
		assertArrayEquals(count, haplotypeModelRandom.getStateFrequencies(), 1e-6);

	}
	
	@Test
	public void testIdentifiable() throws Exception {
		assertNull(haplotypeModel.getId());
		
		haplotypeModel.setId("Test1");
		assertEquals("Test1", haplotypeModel.getId());
		haplotypeModel.setId("Test2");
		assertEquals("Test2", haplotypeModel.getId());
	}
	@Test
	public void testIterable() throws Exception {

		// Iterator()
		Iterator<Taxon> taxonIterator = haplotypeModel.iterator();
		SortedSet<Taxon> attributeSet = new TreeSet<Taxon>();
		while (taxonIterator.hasNext()) {
			Taxon taxon = (Taxon) taxonIterator.next();
			attributeSet.add(taxon);
			taxonIterator.remove();
		}

		SortedSet<Taxon> expectedSet = new TreeSet<Taxon>();
		Collections.addAll(expectedSet, expectedTaxons);

		assertEquals(expectedSet.size(), attributeSet.size());
		for (Iterator<Taxon> iterator = expectedSet.iterator(); iterator.hasNext();) {
			Taxon taxon = iterator.next();
			assertTrue(attributeSet.contains(taxon));
		}
		assertTrue(attributeSet.containsAll(expectedSet));
		
	}
	
	@Test
	public void testAttributable() throws Exception {
		// null test first
		assertNull(haplotypeModel.getAttribute("NULL"));
		assertNull(haplotypeModel.getAttributeNames());
		
		haplotypeModel.setAttribute("A1", 0);
		haplotypeModel.setAttribute("A2", "StringAttribute");
		ArrayList<String> tempList = new ArrayList<String>();
		tempList.add("T1");tempList.add("T2");
		haplotypeModel.setAttribute("A3", tempList);
		int[] tempArray = new int[]{1,2,3,0};
		haplotypeModel.setAttribute("A4", tempArray);
		
		assertEquals(0, haplotypeModel.getAttribute("A1"));
		assertEquals("StringAttribute", haplotypeModel.getAttribute("A2"));
		ArrayList<String> expectedList = new ArrayList<String>();
		expectedList.add("T1");expectedList.add("T2");
		assertEquals(expectedList, haplotypeModel.getAttribute("A3"));
		int[] expectedArray = new int[]{1,2,3,0};
		assertArrayEquals(expectedArray, (int[]) haplotypeModel.getAttribute("A4"));
		assertNull(haplotypeModel.getAttribute("A5"));
		
		Iterator<String> attributeNames = haplotypeModel.getAttributeNames();
		Set<String> attributeSet = new HashSet<String>();
		while (attributeNames.hasNext()) {
			String string = (String) attributeNames.next();
			attributeSet.add(string);
		}

		
		String[] expectedName = new String[]{"A4", "A3", "A1", "A4", "A2"};
		Set<String> expectedSet = new HashSet<String>();
		Collections.addAll(expectedSet, expectedName);

		assertEquals(expectedSet, attributeSet);
			
	}
	
	

	
	
}
