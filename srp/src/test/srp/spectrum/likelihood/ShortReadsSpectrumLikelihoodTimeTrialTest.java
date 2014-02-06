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

//	public static final double ERROR = ShortReadsSpectrumLikelihood.ERROR_RATE;
//	public static final double NOT_ERROR = ShortReadsSpectrumLikelihood.NOT_ERROR_RATE;
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
//		Alignment alignment = DataImporter.importShortReads("/home/sw167/workspaceSrp/snowgoose/srp/unittest/", "5000Srp.fasta");
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

	@Test
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
		String summary = timeTrialOperator(likelihood, op, 1e4);
		System.out.println(summary + "\t" + op.getOperatorName());

	}

//	@Test
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
	
	@Test
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

/*H10

TimeTrial: 5171	0.05171/calculation	DeltaExchangeSingleSpectrumOperator
TimeTrial: 4472	0.4472/calculation	DeltaExchangeMultiSpectrumOperator
TimeTrial: 2539	0.2539/calculation	DeltaExchangeColumnSpectrumOperator
TimeTrial: 505	0.0505/calculation	SwapSingleSpectrumOperator
TimeTrial: 2423	0.2423/calculation	SwapMultiSpectrumOperator
TimeTrial: 3218	0.3218/calculation	DirichletSpectrumOperator
TimeTrial: 644	0.0644/calculation	DirichletAlphaSpectrumOperator
TimeTrial: 7026	7.026/calculation	RecombinationSpectrumOperator
TimeTrial: 2103	0.2103/calculation	RecombineSectionSpectrumOperator

TimeTrial:  	808	8.08/calculation	Full calculation no operator
TimeTrial:	332	0.0332/calculation	Single No store/restore
TimeTrial:	12	1.2E-4/calculation	StoreRestoreOnly


#####################################
## only caluclate delta eachSrpLogLikelihood
TimeTrial: 5815	0.05815/calculation	DeltaExchangeSingleSpectrumOperator
## update eachSrpAt(i) and sum(eachSrp)
TimeTrial: 12199	0.12199/calculation	DeltaExchangeSingleSpectrumOperator


##########################################
## calculate freq*error at each ijk
TimeTrial: 1387	0.1387/calculation	DeltaExchangeSingleSpectrumOperator
TimeTrial: 7323	0.7323/calculation	DeltaExchangeMultiSpectrumOperator
## pre calculate freq*error for each k
TimeTrial: 12199	0.12199/calculation	DeltaExchangeSingleSpectrumOperator
TimeTrial: 5497	0.5497/calculation	DeltaExchangeMultiSpectrumOperator

#########################################
## srpArray vs HashSet on multi, 5000reads
Hashset
TimeTrial: 60128	6.0128/calculation	DeltaExchangeMultiSpectrumOperator
sryArray
TimeTrial: 51010	5.101/calculation	DeltaExchangeMultiSpectrumOperator


##############################
old
TimeTrial:	889	0.0889/calculation	SwapSectionSpectrumOperator //swap index
TimeTrial:	1483	0.1483/calculation	SwapSectionSpectrumOperator // full calculation
*/