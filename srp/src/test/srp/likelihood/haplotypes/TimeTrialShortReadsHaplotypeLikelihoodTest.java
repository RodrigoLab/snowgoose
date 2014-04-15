package test.srp.likelihood.haplotypes;


import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.evolution.OperationType;
import srp.evolution.shortreads.ShortReadMapping;
import srp.evolution.spectrum.Spectrum;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import srp.haplotypes.AlignmentUtils;
import srp.haplotypes.HaplotypeModel;
import srp.likelihood.haplotypes.ShortReadsHaplotypeLikelihood;
import srp.likelihood.spectrum.ShortReadsSpectrumLikelihood;
import srp.likelihood.spectrum.AbstractShortReadsSpectrumLikelihood.DistType;
import srp.operator.haplotypes.AbstractHaplotypeOperator;
import srp.operator.haplotypes.BaseSingleOperator;
import srp.operator.haplotypes.BasesMultiOperator;
import srp.operator.haplotypes.HaplotypeRecombinationOperator;
import srp.operator.haplotypes.HaplotypeSwapSectionOperator;
import srp.operator.spectrum.AbstractSpectrumOperator;
import srp.operator.spectrum.DeltaExchangeColumnSpectrumOperator;
import srp.operator.spectrum.DeltaExchangeMultiSpectrumOperator;
import srp.operator.spectrum.DirichletSpectrumOperator;
import srp.operator.spectrum.RecombinationSpectrumOperator;
import srp.operator.spectrum.RecombineSectionSpectrumOperator;
import srp.operator.spectrum.SwapMultiSpectrumOperator;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.inference.markovchain.MarkovChain;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.math.MathUtils;

/*
 * reflection example
 * 	
ShortReadsHaplotypeLikelihood srpLikelihood = new ShortReadsHaplotypeLikelihood(spectrumModel);
Method method = ShortReadsHaplotypeLikelihood.class.getDeclaredMethod("caluclateStateLogLikelihood", double.class);
method.setAccessible(true);
double log = (Double) method.invoke(srpLikelihood, 0.5);

 */

public class TimeTrialShortReadsHaplotypeLikelihoodTest {

	public static final double ERROR = ShortReadsHaplotypeLikelihood.ERROR_RATE;
	public static final double NOT_ERROR = ShortReadsHaplotypeLikelihood.NOT_ERROR_RATE;
	private static final double THRESHOLD = MarkovChain.EVALUATION_TEST_THRESHOLD;
	
	private static boolean combine = false;
//	private static final double EVALUATION_TEST_THRESHOLD = 1e-8;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	private HaplotypeModel haplotypeModel;
	private ShortReadMapping srpMap;
	private ShortReadsHaplotypeLikelihood likelihood;

	@Before
	public void setUp() throws Exception {

		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "H10_srpGap.fasta");
		ShortReadMapping srpMap = new ShortReadMapping(alignment);
			
//		int spectrumLength = srpMap.getLength();
//		spectrumModel = new SpectrumAlignmentModel(spectrumLength,  6);
//		likelihood = new ShortReadsSpectrumLikelihood(spectrumModel, srpMap, DistType.betaMean);
//
////		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "HaplotypeModelTest_10_srp.fasta");
////		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "H4_srp.fasta");
//		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "SpectrumTest_50srp_200bp.fasta");
//		srpMap = new ShortReadMapping(alignment);
		haplotypeModel = new HaplotypeModel(6, srpMap.getLength());
		likelihood = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void testTimeTrialFull() throws Exception {
		double trial = 100;
		long time1 = System.currentTimeMillis();
		for (int t = 0; t < trial; t++) {
			likelihood.makeDirty();
			likelihood.getLogLikelihood();
		}
		long totalTime = System.currentTimeMillis() - time1;
		System.out.println("TimeTrial: " + totalTime + "\t" + totalTime / trial
				+ "/calculation\tFull calculation no operator");

	}

	@Test
	public void testTimeTrialFullMaster() throws Exception {
		double trial = 100;
		long time1 = System.currentTimeMillis();
		for (int t = 0; t < trial; t++) {
			likelihood.makeDirty();
			likelihood.calculateSrpLikelihoodFullMaster();
		}
		long totalTime = System.currentTimeMillis() - time1;
		System.out.println("TimeTrial: " + totalTime + "\t" + totalTime / trial
				+ "/calculation\tFullMaster calculation no operator");

	}
	
	@Test
	public void testTimeTrialSingle() throws Exception {
		double ite = 5e4;
		for (int i = 0; i < 10; i++) {
			System.out.print("Run "+i +"\t");
		
			AbstractHaplotypeOperator op = new BaseSingleOperator(haplotypeModel);
			String summary = timeTrialOperator(likelihood, op, ite);
			System.out.println(summary );
		}
		System.out.println();
	}
	
	@Test
	public void testTimeTrialMulti() throws Exception {
		double ite = 5e4;
		for (int i = 0; i < 10; i++) {
			System.out.print("Run "+i +"\t");
		
			int bases = 5;
			AbstractHaplotypeOperator op = new BasesMultiOperator(haplotypeModel, bases, CoercionMode.COERCION_OFF);
			String summary = timeTrialOperator(likelihood, op, ite);
			System.out.println(summary +"\t"+ bases);
		}
		System.out.println();
	}	
	
	@Test
	public void testTimeTrialRecombination() throws Exception {
		double ite = 1e4;
		for (int i = 0; i < 10; i++) {
			System.out.print("Run "+i +"\t");
		
			int bases = 5;
//			AbstractHaplotypeOperator op = new HaplotypeRecombinationOperator(haplotypeModel, 0);
			AbstractHaplotypeOperator op = new HaplotypeSwapSectionOperator(haplotypeModel, bases, null);
			String summary = timeTrialOperator(likelihood, op, ite);
			System.out.println(summary +"\t"+ bases);
		}
		System.out.println();
	}
	
	public static String timeTrialOperator(
			ShortReadsHaplotypeLikelihood likelihood2,
			AbstractHaplotypeOperator op, double ite) {
	
		if(combine ){
			return timeTrialOperatorCombine(likelihood2, op, ite);
		}
		else{
			return timeTrialOperatorEach(likelihood2, op, ite);
		}
	}
	
	public static String timeTrialOperatorCombine(
			ShortReadsHaplotypeLikelihood likelihood2,
			AbstractHaplotypeOperator op, double ite) {
		System.gc();
		double scale = 1e3;
		int count = 0;
		long totalTime = 0;
		do {
			try {
				long time1 = System.nanoTime();
//				long time1 = System.nanoTime();
				likelihood2.storeModelState();
				op.doOperation();

				likelihood2.getLogLikelihood();
				double rand = MathUtils.nextDouble();
				if (rand > 0.5) {
					likelihood2.acceptModelState();
				} else {
					likelihood2.restoreModelState();
				}
				totalTime += (System.nanoTime() - time1);
//				totalTime += (System.nanoTime() - time1);

				count++;
			} catch (OperatorFailedException e) {
			}
		} while (count < ite);
		String summary = "TimeTrial: " + totalTime/scale + "\t" + totalTime/scale/ite
				+ "/calculation\t"+ite+" ite.\t"+op.getOperatorName();
		return summary;
	}


	public static String timeTrialOperatorEach(
			ShortReadsHaplotypeLikelihood likelihood2,
			AbstractHaplotypeOperator op, double ite) {
		System.gc();
		double scale = ite*1e3;
		int count = 0;
		long operatorTime = 0;
		long likelihoodTime = 0;
		long storeTime = 0;
		long totalTime = System.nanoTime();

		do {
			try {
				long time0 = System.nanoTime();
				likelihood2.storeModelState();
				storeTime += (System.nanoTime() - time0);
				
				time0 = System.nanoTime();
				op.doOperation();
				operatorTime += (System.nanoTime() - time0);
				
				time0 = System.nanoTime();
				likelihood2.getLogLikelihood();
				likelihoodTime += (System.nanoTime() - time0);
				
				double rand = MathUtils.nextDouble();
				time0 = System.nanoTime();
				if (rand > 0.5) {
					likelihood2.acceptModelState();
				} else {
					likelihood2.restoreModelState();
				}
				storeTime += (System.nanoTime() - time0);
//				totalTime += (System.nanoTime() - time1);

				count++;
			} catch (OperatorFailedException e) {
			}
		} while (count < ite);
		String summary = "TimeTrial: " + operatorTime/scale + "\t" + likelihoodTime/scale
				+ "\t" + storeTime/scale + "\t" + (System.nanoTime()-totalTime)/scale +
				"/calculation\t" + ite + " ite.\t"+op.getOperatorName();

		return summary;
	}

	
}
