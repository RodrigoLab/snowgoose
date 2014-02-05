package test.srp.spectrum.likelihood;


import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.haplotypes.AlignmentMapping;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.likelihood.ShortReadsSpectrumLikelihood;
import srp.spectrum.operator.AbstractSpectrumOperator;
import srp.spectrum.operator.DeltaExchangeColumnSpectrumOperator;
import srp.spectrum.operator.DeltaExchangeMultiSpectrumOperator;
import srp.spectrum.operator.DeltaExchangeSingleSpectrumOperator;
import srp.spectrum.operator.DirichletAlphaSpectrumOperator;
import srp.spectrum.operator.DirichletSpectrumOperator;
import srp.spectrum.operator.RecombinationSpectrumOperator;
import srp.spectrum.operator.RecombineSectionSpectrumOperator;
import srp.spectrum.operator.SwapMultiSpectrumOperator;
import srp.spectrum.operator.SwapSingleSpectrumOperator;
import dr.evolution.alignment.Alignment;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class ShortReadsSpectrumLikelihoodTimeTrialTest {

	public static final double ERROR = ShortReadsSpectrumLikelihood.ERROR_RATE;
	public static final double NOT_ERROR = ShortReadsSpectrumLikelihood.NOT_ERROR_RATE;
	public ShortReadsSpectrumLikelihood likelihood;
	public SpectrumAlignmentModel spectrumModel;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "H10_srp.fasta");
		AlignmentMapping aMap = new AlignmentMapping(alignment);
			
		spectrumModel = new SpectrumAlignmentModel(aMap, 6, 2);
		likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);

	}

	@After
	public void tearDown() throws Exception {
	}


	@Test
	public void testTimeTrialFull() throws Exception {

//		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/ABI/unittest/", "H4_srp.fasta");
//		AlignmentMapping aMap = new AlignmentMapping(alignment);
//			
//		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 4);
//		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double trial = 100;
		long time1 = System.currentTimeMillis();
		for (int t = 0; t < trial; t++) {
			likelihood.makeDirty();
			likelihood.getLogLikelihood();
		}
		long totalTime = System.currentTimeMillis() - time1;
		System.out.println("TimeTrial:  \t"+ totalTime +"\t"+ totalTime/trial +"/calculation\tFull calculation no operator");
		
	}

	@Test
	public void testTimeTrialSingleNoStoreRestore() throws Exception {

//		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/ABI/unittest/", "H4_srp.fasta");
//		AlignmentMapping aMap = new AlignmentMapping(alignment);
//			
//		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 4);
//		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);

		DeltaExchangeSingleSpectrumOperator op = new DeltaExchangeSingleSpectrumOperator(spectrumModel, 0.1, null);

		double trial = 1e4;
		long totalTime = 0;
		int count = 0;
		do{
			try {
				op.doOperation();
				
				long time1 = System.currentTimeMillis();
//				likelihood.makeDirty();
				likelihood.getLogLikelihood();
				totalTime += (System.currentTimeMillis()-time1);
				count++;
				
			} catch (Exception e) {
			}
		}while(count< trial);
		System.out.println("TimeTrial:\t"+ totalTime +"\t"+ totalTime/trial +"/calculation\tSingle No store/restore");
	}


	@Test
	public void testTimeTrialStoreRestoreOnly() throws Exception {

//		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/ABI/unittest/", "H4_srp.fasta");
//		AlignmentMapping aMap = new AlignmentMapping(alignment);
//			
//		SpectrumAlignmentModel spectrumModel = new SpectrumAlignmentModel(aMap, 4);
//		ShortReadsSpectrumLikelihood likelihood = new ShortReadsSpectrumLikelihood(spectrumModel);
		
		double trial = 1e5;
		long totalTime = 0;
		for (int t = 0; t < trial; t++) {
			long time1 = System.currentTimeMillis();
			likelihood.storeModelState();
			likelihood.restoreModelState();
			totalTime += (System.currentTimeMillis()-time1);

		}
		System.out.println("TimeTrial:\t"+ totalTime +"\t"+ totalTime/trial +"/calculation\tStoreRestoreOnly");

	}

	@Test
	public void testTimeTrialAllOperators() throws Exception {
		System.out.println();
		testTimeTrialDeltaSingle();
		testTimeTrialDeltaMulti();
		testTimeTrialDeltaColumn();
		
		testTimeTrialSwapSingle();
		testTimeTrialSwapMulti();
		
		testTimeTrialDirichlet();
		testTimeTrialDirichletAlpha();
		
		testTimeTrialRecombination();
		testTimeTrialRecombineSection();
		System.out.println("");
	}

	
	public void testTimeTrialDeltaSingle() throws Exception {

		AbstractSpectrumOperator op = new DeltaExchangeSingleSpectrumOperator(
				spectrumModel, 0.1, null);
		String summary = timeTrialOperator(likelihood, op, 100000);
		System.out.println(summary + "\t" + op.getOperatorName());

	}

	@Test
	public void testTimeTrialDeltaMulti() throws Exception {

		DeltaExchangeMultiSpectrumOperator op = new DeltaExchangeMultiSpectrumOperator(
				spectrumModel, 0.1, 10, CoercionMode.COERCION_OFF);
		String summary = timeTrialOperator(likelihood, op, 10000);
		System.out.println(summary + "\t" + op.getOperatorName());

	}

	public void testTimeTrialDeltaColumn() throws Exception {

		DeltaExchangeColumnSpectrumOperator op = new DeltaExchangeColumnSpectrumOperator(
				spectrumModel, 0.1, null);
		String summary = timeTrialOperator(likelihood, op, 10000);
		System.out.println(summary + "\t" + op.getOperatorName());

	}

	public void testTimeTrialSwapSingle() throws Exception {

		AbstractSpectrumOperator op = new SwapSingleSpectrumOperator(
				spectrumModel);
		String summary = timeTrialOperator(likelihood, op, 10000);
		System.out.println(summary + "\t" + op.getOperatorName());
	}

	public void testTimeTrialSwapMulti() throws Exception {

		AbstractSpectrumOperator op = new SwapMultiSpectrumOperator(
				spectrumModel, 5, CoercionMode.COERCION_OFF, true);
		String summary = timeTrialOperator(likelihood, op, 10000);
		System.out.println(summary + "\t" + op.getOperatorName());
	}
	
	
	public void testTimeTrialDirichlet() throws Exception {

		AbstractSpectrumOperator op = new DirichletSpectrumOperator(
				spectrumModel, 5, null);
		String summary = timeTrialOperator(likelihood, op, 10000);
		System.out.println(summary + "\t" + op.getOperatorName());
	}

//	@Test
	public void testTimeTrialDirichletAlpha() throws Exception {

		AbstractSpectrumOperator op = new DirichletAlphaSpectrumOperator(
				spectrumModel, 100, null);
		String summary = timeTrialOperator(likelihood, op, 10000);
		System.out.println(summary + "\t" + op.getOperatorName());
	}

	public void testTimeTrialRecombination() throws Exception {

		RecombinationSpectrumOperator op = new RecombinationSpectrumOperator(
				spectrumModel);
		String summary = timeTrialOperator(likelihood, op, 1000);
		System.out.println(summary + "\t" + op.getOperatorName());

	}

	public void testTimeTrialRecombineSection() throws Exception {

		AbstractSpectrumOperator op = new RecombineSectionSpectrumOperator(
				spectrumModel, 10, null);
		String summary = timeTrialOperator(likelihood, op, 10000);
		System.out.println(summary + "\t" + op.getOperatorName());

	}

	public static String timeTrialOperator(
			ShortReadsSpectrumLikelihood likelihood, AbstractSpectrumOperator op, double ite) {
		
		int count = 0;
		long totalTime = 0;
		do {
			try {
				long time1 = System.currentTimeMillis();
				likelihood.storeModelState();
				op.doOperation();

				likelihood.getLogLikelihood();
				double rand = MathUtils.nextDouble();
				if (rand > 0.5) {
					likelihood.acceptModelState();
				} else {
					likelihood.restoreModelState();
				}
				totalTime += (System.currentTimeMillis() - time1);

				count++;
			} catch (OperatorFailedException e) {
			}
		} while (count < ite);
		String summary = "TimeTrial: " + totalTime + "\t" + totalTime / ite
				+ "/calculation";
		return summary;
	}

}

/*
TimeTrial: 1387	0.1387/calculation	DeltaExchangeSingleSpectrumOperator
TimeTrial: 4275	0.4275/calculation	DeltaExchangeMultiSpectrumOperator
TimeTrial: 2941	0.2941/calculation	DeltaExchangeColumnSpectrumOperator
TimeTrial: 1108	0.1108/calculation	SwapSingleSpectrumOperator
TimeTrial: 4377	0.4377/calculation	SwapMultiSpectrumOperator
TimeTrial: 4432	0.4432/calculation	DirichletSpectrumOperator
TimeTrial: 1105	0.1105/calculation	DirichletAlphaSpectrumOperator
TimeTrial: 70782	7.0782/calculation	RecombinationSpectrumOperator
TimeTrial: 2130	0.213/calculation	RecombineSectionSpectrumOperator

TimeTrial:  	6347	63.47/calculation	Full calculation no operator
TimeTrial: 1324	0.1324/calculation	DeltaExchangeSingleSpectrumOperator
TimeTrial:	959	0.0959/calculation	Single No store/restore
TimeTrial:	42	4.2E-4/calculation	StoreRestoreOnly

old
TimeTrial:	889	0.0889/calculation	SwapSectionSpectrumOperator //swap index
TimeTrial:	1483	0.1483/calculation	SwapSectionSpectrumOperator // full calculation
*/