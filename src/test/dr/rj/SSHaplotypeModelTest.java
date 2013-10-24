package test.dr.rj;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotSame;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.haplotypes.AbstractHaplotypeModel;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.Haplotype;
import srp.haplotypes.Operation;
//import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.operator.MultiBasesOperator;
import srp.haplotypes.operator.SingleBaseOperator;
import test.TestUtils;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.math.MathUtils;
import dr.rj.SSHaplotypeModel;


public class SSHaplotypeModelTest {

	private static AlignmentMapping aMap;
	private static Alignment srpAlignment;
	
	private static Taxon[] expectedTaxons;
	private static String[] expectedSequences;
	private static ArrayList<Haplotype> expectedList;
//	private static HashMap<Character, Integer> charToState;

	
	private SSHaplotypeModel haplotypeModel;
	private AbstractHaplotypeModel haplotypeModelRandom;

	
	
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {

		String dir = System.getProperty("user.dir")+File.separatorChar+"unittest"+File.separator;
		srpAlignment = DataImporter.importAlignment(dir, "HaplotypeModelTest_10.fasta");
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
		
//		charToState = new HashMap<Character, Integer>();
//		charToState.put('A', 0);
//		charToState.put('C', 1);
//		charToState.put('G', 2);
//		charToState.put('T', 3);
//		charToState.put('-', 17);
//		charToState.put('*', 17);
//		charToState.put('.', 17);
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		haplotypeModel = new SSHaplotypeModel(aMap, srpAlignment);
		haplotypeModelRandom = new SSHaplotypeModel(aMap, 5);
	}

	@After
	public void tearDown() throws Exception {
	}
	
	
	@Test
	public void testInit() throws Exception {
		
		assertEquals(16, haplotypeModel.getHaplotypeCount());
		assertEquals(16, haplotypeModel.getSequenceCount());
		assertEquals(16, haplotypeModel.getTaxonCount());
		assertEquals(100, haplotypeModel.getHaplotypeLength());
		for (int i = 0; i < haplotypeModel.getTrueHaplotypeCount(); i++) {
			Taxon expected = new Taxon("r"+(i+1)+".1");
			assertEquals(expected, haplotypeModel.getTaxon(i));
			System.out.println(i +"\t"+ haplotypeModel.getHaplotypeString(i));
		}
		for (int i = haplotypeModel.getTrueHaplotypeCount(); i < haplotypeModel
				.getHaplotypeCount(); i++) {
			Taxon expected = new Taxon(SSHaplotypeModel.TAXON_PREFIX+i);
			assertEquals(expected, haplotypeModel.getTaxon(i));
			System.err.println(i +"\t"+ haplotypeModel.getHaplotypeString(i));
		}

		assertEquals(16, haplotypeModelRandom.getHaplotypeCount());
		assertEquals(16, haplotypeModelRandom.getSequenceCount());
		assertEquals(16, haplotypeModelRandom.getTaxonCount());
		assertEquals(100, haplotypeModelRandom.getHaplotypeLength());
		for (int i = 0; i < haplotypeModelRandom.getHaplotypeCount(); i++) {
			Taxon expected = new Taxon(SSHaplotypeModel.TAXON_PREFIX+i);
			assertEquals(expected, haplotypeModelRandom.getTaxon(i));
			System.out.println(i +"\t"+ haplotypeModelRandom.getHaplotypeString(i));
		}
	}
	
	@Test
	public void testGetAlignment() throws Exception {
		
	
		String[] newHaplotypes = expectedSequences;

		for (int i = 0; i < newHaplotypes.length; i++) {
			assertEquals(haplotypeModel.getHaplotypeString(i), haplotypeModel.getAlignedSequenceString(i)) ;
			assertEquals(newHaplotypes[i], haplotypeModel.getHaplotypeString(i));
			for (int j = 0; j < newHaplotypes.length; j++) {
				char expected = newHaplotypes[i].charAt(j); 
				assertEquals(expected,haplotypeModel.getHaplotypeCharAt(i,j));
			}
			
		}
	}
	
	@Test
	public void testReject() throws Exception {
		
//		haplotypeModel = new SSHaplotypeModel(aMap, 12);		
		int trueHaplotypeCount = haplotypeModel.getTrueHaplotypeCount();
		char[] temp = new char[haplotypeModel.getHaplotypeLength()];
		Arrays.fill(temp, DataType.GAP_CHARACTER);
		String tempSeq = String.valueOf(temp);
		
		
		AlignmentMapping aMap = haplotypeModel.getAlignmentMapping(); 
		for (int i = 0; i < 1000; i++) {

			int[] posChar = aMap.getNextBase();
			
			haplotypeModel.swapHaplotypeSingleBase(Operation.SWAPSINGLE, posChar);
			haplotypeModel.reject();
		}
		for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
			if(i<trueHaplotypeCount){
				assertEquals(haplotypeModel.getAlignedSequenceString(i), srpAlignment.getAlignedSequenceString(i));
			}
			else{
				assertEquals(haplotypeModel.getAlignedSequenceString(i), tempSeq);
			}
		}

		int[][] allPosChars = new int[2][haplotypeModel.getHaplotypeLength()];
		for (int t = 0; t < 1000; t++) {
			
			for (int i = 0; i < allPosChars.length; i++) {
				Arrays.fill(allPosChars[i], -1);
			}
			int swapLength = MathUtils.nextInt(haplotypeModel.getHaplotypeLength());
			for (int i = 0; i < swapLength; i++) {
				int[] posChar = aMap.getNextBase();
				allPosChars[0][posChar[0]] = posChar[1];
			}
			haplotypeModel.swapHaplotypeMultiBases(Operation.SWAPMULTI, allPosChars);

			
			haplotypeModel.reject();
		}
		for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
			if(i<trueHaplotypeCount){
				assertEquals(haplotypeModel.getAlignedSequenceString(i), srpAlignment.getAlignedSequenceString(i));
			}
			else{
				assertEquals(haplotypeModel.getAlignedSequenceString(i), tempSeq);
			}
		}
	
	}
	

	@Test
	public void testSetupAlignment() throws Exception {

		haplotypeModel = new SSHaplotypeModel(aMap, 5);
		for (int i = 0; i < haplotypeModel.getTrueHaplotypeCount(); i++) {
			assertNotEquals(haplotypeModel.getHaplotypeString(i), srpAlignment.getAlignedSequenceString(i));
			assertNotSame(haplotypeModel.getHaplotypeString(i), srpAlignment.getSequence(i));
			assertNotEquals(haplotypeModel.getHaplotypeString(i), expectedSequences[0] );
			assertNotSame(haplotypeModel.getHaplotypeString(i), expectedSequences[0]);
		
			assertNotSame(haplotypeModel.getHaplotypeString(i), srpAlignment.getSequence(i));
			assertNotEquals(haplotypeModel.getHaplotypeString(i), expectedSequences[0] );
			assertNotSame(haplotypeModel.getHaplotypeString(i), expectedSequences[0]);
		}	
	}

	@Test
	public void testHaplotypeGetNextBaseFrequency() throws Exception {


		Parameter frequency = new Parameter.Default(new double[]{0.1,0.25,0.60,0.05});

		int[] count = new int['T'+1]; 
		int ite = 1000;
    	for (int i = 0; i < ite; i++) {
    		int[] posChar = haplotypeModel.getNextBaseFrequency(frequency);
    		count[posChar[1]]++;
		}
    	TestUtils.assertExpectationRange(count['A'], frequency.getParameterValue(0)*ite, 0.05*ite);
    	TestUtils.assertExpectationRange(count['C'], frequency.getParameterValue(1)*ite, 0.05*ite);
    	TestUtils.assertExpectationRange(count['G'], frequency.getParameterValue(2)*ite, 0.05*ite);
    	TestUtils.assertExpectationRange(count['T'], frequency.getParameterValue(3)*ite, 0.05*ite);

    	
    	
    	double[] logFreq = new double[4];
    	for (int i = 0; i < logFreq.length; i++) {
			logFreq[i] = Math.log(frequency.getParameterValue(i));
		}
    	for (int i = 0; i < logFreq.length; i++) {
			for (int j = 0; j < logFreq.length; j++) {
				double expected = logFreq[i]-logFreq[j];
				assertEquals(expected, haplotypeModel.getLogqFrequencyStates(i,j), 0);
			}
		}
    	
    	frequency = new Parameter.Default(new double[]{0.1, 0.2, 0.3, 0.4});
    	haplotypeModel.getNextBaseFrequency(frequency);
     	
    	for (int i = 0; i < logFreq.length; i++) {
			logFreq[i] = Math.log(frequency.getParameterValue(i));
		}
    	for (int i = 0; i < logFreq.length; i++) {
			for (int j = 0; j < logFreq.length; j++) {
				double expected = logFreq[i]-logFreq[j];
				assertEquals(expected, haplotypeModel.getLogqFrequencyStates(i,j), 0);
			}
		}
    	count = new int['T'+1]; 
    	for (int i = 0; i < ite; i++) {
    		int[] posChar = haplotypeModel.getNextBaseFrequency(frequency);
    		count[posChar[1]]++;
		}
    	TestUtils.assertExpectationRange(count['A'], frequency.getParameterValue(0)*ite, 0.05*ite);
    	TestUtils.assertExpectationRange(count['C'], frequency.getParameterValue(1)*ite, 0.05*ite);
    	TestUtils.assertExpectationRange(count['G'], frequency.getParameterValue(2)*ite, 0.05*ite);
    	TestUtils.assertExpectationRange(count['T'], frequency.getParameterValue(3)*ite, 0.05*ite);

    	
    	
	}
	
	
	
	
}
