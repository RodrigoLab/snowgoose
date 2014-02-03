package test.srp.spectrum;


import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;

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

import org.hamcrest.CoreMatchers;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import srp.core.DataImporter;
import srp.dr.evolution.datatype.ShortReads;
import srp.haplotypes.AlignmentMapping;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
//import srp.haplotypes.Haplotype;
//import srp.haplotypes.HaplotypeModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

public class AbstractSpectrumAlignmentModelTest {


	private static AlignmentMapping aMap;
	private static Alignment srpAlignment;
	
	private static Taxon[] expectedTaxons;
	private static String[] expectedSequences;
	private static ArrayList<Spectrum> expectedList;
	private static HashMap<Character, Integer> charToState;
	
	private SpectrumAlignmentModel spectrumModel;
	private SpectrumAlignmentModel spectrumModelRandom;

	
	
	
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
//			char[][] expectedMatrix = new char[matrixS.length][matrixS[0].length()];
		expectedList = new ArrayList<Spectrum>();
		for (int i = 0; i < expectedSequences.length; i++) {
			expectedTaxons[i] = new Taxon("taxa_"+i);
//			Spectrum h = new Spectrum(new Sequence(expectedSequences[i]));
//			expectedList.add(h);
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
		spectrumModel = new SpectrumAlignmentModel(aMap, srpAlignment);
		spectrumModelRandom = new SpectrumAlignmentModel(aMap, 5);
	}

	@After
	public void tearDown() throws Exception {
	}


	@Test
	public void testHaplotypes() throws Exception {
		
		for (int i = 0; i < expectedList.size(); i++) {
		
//			assertEquals(expectedSequences[i], spectrumModel.getAlignedSequenceString(i));
//			assertEquals(expectedSequences[i], spectrumModel.getUnalignedSequenceString(i));
//			assertEquals(expectedSequences[i], haplotypeModel.getHaplotype(i).getSequenceString());
			assertEquals(expectedTaxons[i], spectrumModel.getSpectrum(i).getTaxon());
			
//			assertEquals(haplotypeModel.getSequence(i), haplotypeModel.getHaplotype(i));
//				expectedMatrix[i] = matrixS[i].toCharArray();
			for (int j = 0; j < expectedSequences[i].length(); j++) {
				char expected = expectedSequences[i].charAt(j); 
//				assertEquals(expected,haplotypeModel.getHaplotypeCharAt(i,j));
			}
		}
//			assertArrayEquals(expectedMatrix, haplotypeModel.getCharMatrix());
		
	}
	
	
	@Test
	public void testSequenceList() throws Exception {
		// getSequence(int)
		// in testHaplotypes()
		
		// getSequenceCount()
		// in testInit()
		
		// getSequenceAttribute(int, String)
		// setSequenceAttribute(int, String, Object)
//		for (int i = 0; i < spectrumModel.getSpectrumCount(); i++) {
//			assertNull(spectrumModel.getSequenceAttribute(i, "NULL"));
//			spectrumModel.setSequenceAttribute(i, "NAME", expectedTaxons[i].getId());
//		}
//		for (int i = 0; i < spectrumModel.getSpecturmCount(); i++) {
//			assertEquals(("r"+(i+1)+".1"), spectrumModel.getSequenceAttribute(i, "NAME"));
//			assertNull(spectrumModel.getSequenceAttribute(i, "NULL"));
//		}
		
	}
	@Test
	public void testTaxonList() throws Exception {

		// getTaxonCount()
		assertEquals(10, spectrumModel.getTaxonCount());
		assertEquals(5, spectrumModelRandom.getTaxonCount());

		// getTaxon(int)
		// getTaxonId(int)
		// getTaxonIndex(String)
		// getTaxonIndex(Taxon)
		
//		for (int i = 0; i < expectedTaxons.length; i++) {
//			String expectedID = expectedTaxons[i].getId();
//			assertEquals(expectedTaxons[i], spectrumModel.getTaxon(i));
//			assertEquals(expectedID, spectrumModel.getTaxonId(i));
//			assertEquals(i, spectrumModel.getTaxonIndex(expectedID));
//			assertEquals(i, spectrumModel.getTaxonIndex(spectrumModel.getTaxon(i)));
//			
//		}
//		assertEquals(-1, spectrumModel.getTaxonIndex("EMPTY"));
		assertEquals(-1, spectrumModel.getTaxonIndex(new Taxon("E")));
		assertEquals(-1, spectrumModel.getTaxonIndex(expectedTaxons[0]));
		
		List<Taxon> expectedList = new ArrayList<Taxon>();
		for (int i = 0; i < spectrumModelRandom.getTaxonCount(); i++) {
			String expectedID = SpectrumAlignmentModel.TAXON_PREFIX+i;
			Taxon expectedTaxon = new Taxon(expectedID);
			expectedList.add(expectedTaxon);
			assertEquals(expectedTaxon, spectrumModelRandom.getTaxon(i));
			assertEquals(expectedID, spectrumModelRandom.getTaxonId(i));
			assertEquals(i, spectrumModelRandom.getTaxonIndex(expectedID));
			assertEquals(i, spectrumModelRandom.getTaxonIndex(spectrumModelRandom.getTaxon(i)));
		}
		assertEquals(-1, spectrumModelRandom.getTaxonIndex("EMPTY"));
		assertEquals(-1, spectrumModelRandom.getTaxonIndex(expectedTaxons[0]));
		
		
		// asList()		
//		assertEquals(Arrays.asList(expectedTaxons), spectrumModel.asList());
		assertEquals(expectedList, spectrumModelRandom.asList());

		
		// getTaxonAttribute(int, String)
		for (int i = 0; i < spectrumModelRandom.getTaxonCount(); i++) {
			Taxon taxon = spectrumModelRandom.getTaxon(i);
			taxon.setAttribute("NAME", expectedTaxons[i].getId());
			taxon.setAttribute("INDEX", i);
		}
		for (int i = 0; i < spectrumModelRandom.getTaxonCount(); i++) {
			assertEquals(("taxa_"+i), spectrumModelRandom.getTaxonAttribute(i, "NAME"));
			assertEquals(i, spectrumModelRandom.getTaxonAttribute(i, "INDEX"));
			assertNull(spectrumModelRandom.getTaxonAttribute(i, "NULL"));
			
		}
		
		SimpleAlignment alignmentNoTaxon = new SimpleAlignment();
		alignmentNoTaxon.setDataType(ShortReads.INSTANCE);
		Sequence h;
		for (int i = 0; i < 5; i++) {
//			if (i == 0) {
//				h = AbstractSpectrumAlignmentModelTest.expectedList.get(0);
//			} else {
				h = new Sequence(expectedSequences[i]);
//			}/

			alignmentNoTaxon.addSequence(h);

		}
		
		spectrumModelRandom = new SpectrumAlignmentModel(aMap, alignmentNoTaxon);
		assertEquals(5, spectrumModelRandom.getTaxonCount());
		for (int i = 0; i < spectrumModelRandom.getTaxonCount(); i++) {
			
			assertEquals("taxa_"+i,spectrumModelRandom.getTaxon(i).getId());
			
			assertNull(spectrumModelRandom.getTaxonAttribute(i, "INDEX"));
//			spectrumModelRandom.setSequenceAttribute(i, "INDEX", i+i);
			
//			assertEquals(i+i, spectrumModelRandom.getSequenceAttribute(i,"INDEX"));
//			assertEquals(i, spectrumModelRandom.getTaxonAttribute(i, "INDEX"));
			
		}
		
		for (int i = 0; i < spectrumModelRandom.getTaxonCount(); i++) {

			try {
				spectrumModelRandom.getTaxonId(i);
				
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
			assertEquals(expected, spectrumModel.getPatternIndex(expected));
			assertEquals(expected, spectrumModelRandom.getPatternIndex(expected));
		}
		
		// getSiteCount()
		assertEquals(100, spectrumModel.getSiteCount());
		assertEquals(100, spectrumModelRandom.getSiteCount());
		
		// getSitePattern(int)
		int[] expectedArray = new int[10];
		for (int site = 0; site < expectedSequences[0].length(); site++) {
			for (int i = 0; i < expectedArray.length; i++) {
				expectedArray[i] = charToState.get(expectedSequences[i].charAt(site));
			}
//			assertArrayEquals(expectedArray, spectrumModel.getSitePattern(site));
		}
		Arrays.fill(expectedArray, 17);
//		assertArrayEquals(expectedArray, spectrumModel.getSitePattern(5));
//		assertArrayEquals(expectedArray, spectrumModel.getSitePattern(101));
		
		

		
		// getState(int, int)
//		for (int h = 0; h < haplotypeModelRandom.getHaplotypeCount(); h++) {
//			String hap = haplotypeModelRandom.getHaplotypeString(h);
//			assertEquals(hap.length(), haplotypeModelRandom.getHaplotypeLength());
//			assertEquals(-1, hap.indexOf('.'));
//			for (int i = 0; i < hap.length(); i++) {
//				int expected = charToState.get(hap.charAt(i));
//				assertEquals(expected, haplotypeModelRandom.getState(h, i));
//				assertEquals(expected, haplotypeModelRandom.getPatternState(h, i));
//			}
//			int expected = 17; 
//			assertEquals(expected, haplotypeModelRandom.getState(h, hap.length()));
//			assertEquals(expected, haplotypeModelRandom.getPatternState(h, hap.length()));
//		}
		
//		for (int h = 0; h < haplotypeModel.getHaplotypeCount(); h++) {
//			String hap = expectedSequences[h];
//			assertEquals(hap.length(), haplotypeModel.getHaplotypeLength());
//			for (int i = 0; i < hap.length(); i++) {
//				int expected = charToState.get(hap.charAt(i));
//				assertEquals(expected, haplotypeModel.getState(h, i));
//				assertEquals(expected, haplotypeModel.getPatternState(h, i));
//			}
//			int expected = 17;
//			assertEquals(expected, haplotypeModel.getState(h, hap.length()));
//			assertEquals(expected, haplotypeModel.getPatternState(h, hap.length()));
//		}

	}

	@Test
	public void testPatternList() throws Exception {
		// getDataType()
		assertEquals(Nucleotides.INSTANCE, spectrumModel.getDataType());
		assertEquals(Nucleotides.class, spectrumModel.getDataType().getClass());
		assertEquals(Nucleotides.DESCRIPTION, spectrumModel.getDataType().getDescription());

		// getPattern(int)
		// in testSiteList()
		
		// getPatternCount()
		assertEquals(100, spectrumModel.getPatternCount());
		assertEquals(100, spectrumModelRandom.getPatternCount());
		
		// getPatternLength()
		assertEquals(10, spectrumModel.getPatternLength());
		assertEquals(5, spectrumModelRandom.getPatternLength());
		
		// getPatternState(int, int)
		// in testSiteList()
		
		
		// getPatternWeight(int)
		for (int i = 0; i < 100; i++) {
			assertEquals(1, spectrumModel.getPatternWeight(i), 0);
			assertEquals(1, spectrumModelRandom.getPatternWeight(i), 0);
		}
		
		// getPatternWeights()
		double[] expectedArray = new double[spectrumModel.getPatternCount()];
		Arrays.fill(expectedArray, 1.0);
		assertArrayEquals(expectedArray, spectrumModel.getPatternWeights(), 0);
		assertArrayEquals(expectedArray, spectrumModelRandom.getPatternWeights(), 0);
		
		// getStateCount()
		assertEquals(4, spectrumModel.getStateCount());
		assertEquals(4, spectrumModelRandom.getStateCount());
		
	}
	
	@Test
	public void testIdentifiable() throws Exception {
		assertNull(spectrumModel.getId());
		
		spectrumModel.setId("Test1");
		assertEquals("Test1", spectrumModel.getId());
		spectrumModel.setId("Test2");
		assertEquals("Test2", spectrumModel.getId());
	}
	@Test
	public void testIterable() throws Exception {

		// Iterator()
		Iterator<Taxon> taxonIterator = spectrumModel.iterator();
		SortedSet<Taxon> attributeSet = new TreeSet<Taxon>();
		while (taxonIterator.hasNext()) {
			Taxon taxon = (Taxon) taxonIterator.next();
			attributeSet.add(taxon);
			System.out.println(taxon);
			taxonIterator.remove();
		}

		SortedSet<Taxon> expectedSet = new TreeSet<Taxon>();
		Collections.addAll(expectedSet, expectedTaxons);

		assertEquals(expectedSet.size(), attributeSet.size());
		for (Iterator<Taxon> iterator = expectedSet.iterator(); iterator.hasNext();) {
			Taxon taxon = iterator.next();
			System.out.println(taxon);
			assertTrue(attributeSet.contains(taxon));
		}
		assertTrue(attributeSet.containsAll(expectedSet));
		
	}
	
	@Test
	public void testAttributable() throws Exception {
		// null test first
		assertNull(spectrumModel.getAttribute("NULL"));
		assertNull(spectrumModel.getAttributeNames());
		
		spectrumModel.setAttribute("A1", 0);
		spectrumModel.setAttribute("A2", "StringAttribute");
		ArrayList<String> tempList = new ArrayList<String>();
		tempList.add("T1");tempList.add("T2");
		spectrumModel.setAttribute("A3", tempList);
		int[] tempArray = new int[]{1,2,3,0};
		spectrumModel.setAttribute("A4", tempArray);
		
		assertEquals(0, spectrumModel.getAttribute("A1"));
		assertEquals("StringAttribute", spectrumModel.getAttribute("A2"));
		ArrayList<String> expectedList = new ArrayList<String>();
		expectedList.add("T1");expectedList.add("T2");
		assertEquals(expectedList, spectrumModel.getAttribute("A3"));
		int[] expectedArray = new int[]{1,2,3,0};
		assertArrayEquals(expectedArray, (int[]) spectrumModel.getAttribute("A4"));
		assertNull(spectrumModel.getAttribute("A5"));
		
		Iterator<String> attributeNames = spectrumModel.getAttributeNames();
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
