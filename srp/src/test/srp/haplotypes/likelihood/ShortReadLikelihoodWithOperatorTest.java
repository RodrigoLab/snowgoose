package test.srp.haplotypes.likelihood;

import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.HaplotypeModelUtils;
import srp.haplotypes.Operation;
import srp.haplotypes.likelihood.ShortReadLikelihood;
import srp.haplotypes.operator.AbstractBasesMultiOperator;
import srp.haplotypes.operator.AbstractBaseSingleOperator;
import srp.haplotypes.operator.BaseSingleEmpiricalOperator;
import srp.haplotypes.operator.BaseSingleFrequencyOperator;
import srp.haplotypes.operator.BaseSingleOperator;
import srp.haplotypes.operator.BaseSingleUniformOperator;
import srp.haplotypes.operator.BasesMultiEmpiricalOperator;
import srp.haplotypes.operator.BasesMultiOperator;
import srp.haplotypes.operator.BasesMultiUniformOperator;
import srp.haplotypes.operator.ColumnOperator;
import srp.haplotypes.operator.HaplotypeRecombinationOperator;
import srp.haplotypes.operator.HaplotypeSwapSectionOperator;
import dr.evolution.alignment.Alignment;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.math.MathUtils;

public class ShortReadLikelihoodWithOperatorTest {

	
	private HaplotypeModel haplotypeModelH4;
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {

	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {


		String dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/";
	
		String trueAlignmentFile = "H4_haplotypes.phyml";
		String shortReadFile = "H4_srp.fasta";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		
		Alignment shortReads = dataImporter.importShortReads(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
		
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		haplotypeModelH4 = new HaplotypeModel(aMap, trueAlignment.getSequenceCount());
	
	}

	@After
	public void tearDown() throws Exception {
	}

	
	@Test
	public void testCalculateLikelihoodColumn() throws Exception {
		
		Parameter freqs = new Parameter.Default("frequency", haplotypeModelH4.getStateFrequencies());
		MCMCOperator op = new ColumnOperator(haplotypeModelH4, haplotypeModelH4.getHaplotypeCount(), freqs, null);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, ColumnOperator.OP);
	}
	@Test
	public void testCalculateLikelihoodSingleBase() throws Exception {
		
		MCMCOperator op = new BaseSingleOperator(haplotypeModelH4, 0);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, AbstractBaseSingleOperator.OP);
		
	}
	@Test
	public void testCalculateLikelihoodSingleBaseUniform() throws Exception {
		
		MCMCOperator op = new BaseSingleUniformOperator(haplotypeModelH4, 0);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, AbstractBaseSingleOperator.OP);
	}
	@Test
	public void testCalculateLikelihoodSingleBaseEmpirical() throws Exception {
		
		MCMCOperator op = new BaseSingleEmpiricalOperator(haplotypeModelH4, 0);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, AbstractBaseSingleOperator.OP);
	}
	
	@Test
	public void testCalculateLikelihoodSingleBaseFrequency() throws Exception {
			

		Parameter freqs = new Parameter.Default("frequency", haplotypeModelH4.getStateFrequencies());
		MCMCOperator op = new BaseSingleFrequencyOperator(haplotypeModelH4, freqs);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, AbstractBaseSingleOperator.OP);
		
		
	}
	@Test
	public void testCalculateLikelihoodMultiBases() throws Exception {
		
		MCMCOperator op = new BasesMultiOperator(haplotypeModelH4, 50, CoercionMode.COERCION_OFF);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, AbstractBasesMultiOperator.OP);
	}
	@Test
	public void testCalculateLikelihoodMultiBasesEmpirical() throws Exception {
		
		MCMCOperator op = new BasesMultiEmpiricalOperator(haplotypeModelH4, 50, CoercionMode.COERCION_OFF);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, AbstractBasesMultiOperator.OP);
	}
	@Test
	public void testCalculateLikelihoodMultiBasesUniform() throws Exception {
		
		MCMCOperator op = new BasesMultiUniformOperator(haplotypeModelH4, 50, CoercionMode.COERCION_OFF);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, AbstractBasesMultiOperator.OP);
	}
	@Test
	public void testCalculateLikelihoodSwapSectionRecombination() throws Exception {
		
		MCMCOperator op = new HaplotypeRecombinationOperator(haplotypeModelH4, 0);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, HaplotypeRecombinationOperator.OP);
	}
	@Test
	public void testCalculateLikelihoodSwapSection() throws Exception {
		
		MCMCOperator op = new HaplotypeSwapSectionOperator(haplotypeModelH4, 50, null);
		runTestCalculateSrpLikelihoodOperators(haplotypeModelH4, op, HaplotypeSwapSectionOperator.OP);
	}
	
	public static void runTestCalculateSrpLikelihoodOperators(HaplotypeModel haplotypeModel, MCMCOperator op, Object expectedOperation) throws Exception {
	
		
		ShortReadLikelihood srL = new ShortReadLikelihood(haplotypeModel);
		assertEquals(Operation.NONE, srL.getOperation());
		assertEquals(Operation.NONE, haplotypeModel.getSwapInfo().getOperation());
		for (int i = 0; i < 100; i++) {
	    	srL.storeModelState();
	        op.operate();
	
	        double score = srL.getLogLikelihood();
	        assertEquals(expectedOperation, srL.getOperation());
	        assertEquals(expectedOperation, haplotypeModel.getSwapInfo().getOperation());
	        
	        HaplotypeModel duplicateHaplotypeModel = HaplotypeModelUtils.copyHaplotypeModel(haplotypeModel);
	        ShortReadLikelihood srLFull = new ShortReadLikelihood(duplicateHaplotypeModel);
	        double expected = srLFull.getLogLikelihood();
	        
	        assertEquals(Operation.NONE, srLFull.getOperation());
	        assertEquals(Operation.NONE, duplicateHaplotypeModel.getSwapInfo().getOperation());
	        assertEquals(expected, score, 0);
	
	        boolean accept = MathUtils.nextBoolean();    		
	        if (accept) {
	            op.accept(0);
	            srL.acceptModelState();
	        } else {
	            op.reject();
	            srL.restoreModelState();
	        }
	    }
		
		
	}
}
