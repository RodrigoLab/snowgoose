package test.srp.spectrum.likelihood;


import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import srp.spectrum.SpectrumOperationRecord;
import srp.spectrum.likelihood.ShortReadsSpectrumLikelihood;
import srp.spectrum.operator.AbstractSpectrumOperator;
import srp.spectrum.operator.DeltaExchangeColumnSpectrumOperator;
import srp.spectrum.operator.DeltaExchangeMultiSpectrumOperator;
import srp.spectrum.operator.RecombinationSpectrumOperator;
import srp.spectrum.operator.DeltaExchangeSingleSpectrumOperator;
import srp.spectrum.operator.RecombineSectionSpectrumOperator;
import srp.spectrum.operator.SwapSingleSpectrumOperator;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class ShortReadsSpectrumLikelihoodTimeTrialTest {

	public static final double ERROR = ShortReadsSpectrumLikelihood.ERROR_RATE;
	public static final double NOT_ERROR = ShortReadsSpectrumLikelihood.NOT_ERROR_RATE;
	private ShortReadsSpectrumLikelihood likelihood;
	private SpectrumAlignmentModel spectrumModel;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/ABI/unittest/", "H4_srp.fasta");
		AlignmentMapping aMap = new AlignmentMapping(alignment);
			
		spectrumModel = new SpectrumAlignmentModel(aMap, 4);
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
	public void testTimeTrialDeltaSingle() throws Exception {

		AbstractSpectrumOperator op = new DeltaExchangeSingleSpectrumOperator(
				spectrumModel, 0.1, null);
		String summary = timeTrialOperator(likelihood, op);
		System.out.println(summary + "\t" + op.getOperatorName());

	}
	

	@Test
	public void testTimeTrialSwapSingle() throws Exception {

		AbstractSpectrumOperator op = new SwapSingleSpectrumOperator(
				spectrumModel, null);
		String summary = timeTrialOperator(likelihood, op);
		System.out.println(summary + "\t" + op.getOperatorName());

	}

	@Test
	public void testTimeTrialMulti() throws Exception {

		DeltaExchangeMultiSpectrumOperator op = new DeltaExchangeMultiSpectrumOperator(
				spectrumModel, 0.1, 5, null);
		String summary = timeTrialOperator(likelihood, op);
		System.out.println(summary + "\t" + op.getOperatorName());

	}

	@Test
	public void testTimeTrialColumn() throws Exception {

		DeltaExchangeColumnSpectrumOperator op = new DeltaExchangeColumnSpectrumOperator(
				spectrumModel, 0.1, null);
		String summary = timeTrialOperator(likelihood, op);
		System.out.println(summary + "\t" + op.getOperatorName());

	}

	@Test
	public void testTimeTrialRecombination() throws Exception {

		RecombinationSpectrumOperator op = new RecombinationSpectrumOperator(
				spectrumModel, null);
		String summary = timeTrialOperator(likelihood, op);
		System.out.println(summary + "\t" + op.getOperatorName());

	}
	@Test
	public void testTimeTrialSwapSection() throws Exception {

		AbstractSpectrumOperator op = new RecombineSectionSpectrumOperator(
				spectrumModel, 10, null);
		String summary = timeTrialOperator(likelihood, op);
		System.out.println(summary + "\t" + op.getOperatorName());

	}

	private static String timeTrialOperator(ShortReadsSpectrumLikelihood likelihood,
			AbstractSpectrumOperator op) {
		double trial = 1e4;
		int count = 0;
		long totalTime = 0;
		do{
			try {
//				op.doOperation();
				
				long time1 = System.currentTimeMillis();
				likelihood.storeModelState();
//				likelihood.makeDirty();
				op.doOperation();
				

				likelihood.getLogLikelihood();
				double rand = MathUtils.nextDouble();
				if(rand>0.5){
					likelihood.acceptModelState();
				}
				else{
					likelihood.restoreModelState();
				}
				totalTime += (System.currentTimeMillis()-time1);

				count++;
			} catch (OperatorFailedException e) {
			}
		}while(count< trial);
		String summary = "TimeTrial:\t"+ totalTime +"\t"+ totalTime/trial +"/calculation";
		return summary;
	}

}


/*
TimeTrial:	1006	0.1006/calculation	ColumnSpectrumDeltaExchangeOperator
TimeTrial:	393	0.0393/calculation	SingleSpectrumDeltaExchangeOperator
TimeTrial:  	1858	18.58/calculation	Full calculation no operator
TimeTrial:	323	0.0323/calculation	Single No store/restore
TimeTrial:	678	0.0678/calculation	SwapSectionSpectrumOperator
TimeTrial:	29	2.9E-4/calculation	StoreRestoreOnly
TimeTrial:	20315	2.0315/calculation	RecombinationSpectrumOperator
TimeTrial:	1552	0.1552/calculation	MultiSpectrumDeltaExchangeOperator


TimeTrial:	889	0.0889/calculation	SwapSectionSpectrumOperator //swap index
TimeTrial:	1483	0.1483/calculation	SwapSectionSpectrumOperator // full calculation
*/