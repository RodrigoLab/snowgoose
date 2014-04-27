package test.srp.likelihood.haplotypes;


import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.evolution.haplotypes.HaplotypeModel;
import srp.evolution.shortreads.ShortReadMapping;
import srp.likelihood.haplotypes.ShortReadsHaplotypeLikelihood;
import srp.operator.haplotypes.AbstractHaplotypeOperator;
import srp.operator.haplotypes.BaseSingleOperator;
import srp.operator.haplotypes.BasesMultiOperator;
import srp.operator.haplotypes.ColumnOperator;
import srp.operator.haplotypes.HaplotypeSwapSectionOperator;
import dr.evolution.alignment.Alignment;
import dr.inference.markovchain.MarkovChain;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
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
		int noRun = 10;
		
		AbstractHaplotypeOperator op = new BaseSingleOperator(haplotypeModel);
		timeTrialOperator(likelihood, op, ite, noRun);

	}
	
	
	@Test
	public void testTimeTrialColumn() throws Exception {
		double ite = 5e4;
		int noRun = 10;
		
		AbstractHaplotypeOperator op = new ColumnOperator(haplotypeModel, haplotypeModel.getHaplotypeCount(), null);
		String summary = timeTrialOperator(likelihood, op, ite, noRun);
		
	}
	
	@Test
	public void testTimeTrialMulti() throws Exception {
		double ite = 5e4;
		int run = 10;
		
		int bases = 5;
		AbstractHaplotypeOperator op = new BasesMultiOperator(haplotypeModel, bases, CoercionMode.COERCION_OFF);
		
		timeTrialOperator(likelihood, op, ite, run, bases);
		
	}	
	
	@Test
	public void testTimeTrialRecombination() throws Exception {
		double ite = 5e4;
		int run = 10;

		int bases = 5;
		AbstractHaplotypeOperator op = new HaplotypeSwapSectionOperator(haplotypeModel, bases, null);
		
		timeTrialOperator(likelihood, op, ite, run, bases );
		
	}
	
	public static String timeTrialOperator(
			ShortReadsHaplotypeLikelihood likelihood2,
			AbstractHaplotypeOperator op, double ite, int runs, int... options) {
	
		String summary = null;
		for (int i = 0; i < runs; i++) {
			System.out.print("Run "+i +"\t");
			if(combine ){
				summary = timeTrialOperatorCombine(likelihood2, op, ite);
			}
			else{
				summary = timeTrialOperatorEach(likelihood2, op, ite);
			}
			System.out.println(summary +"\t"+ Arrays.toString(options));
		}
		System.out.println();
		return summary;
		
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

/*
Run 0	TimeTrial: 0.66236218	29.22460532	7.71846686	38.202374/calculation	50000.0 ite.	ColumnOperator	[]
Run 1	TimeTrial: 0.41234174	29.10338322	6.5255292	36.2145474/calculation	50000.0 ite.	ColumnOperator	[]
Run 2	TimeTrial: 0.47637204	32.22797204	7.02296978	39.91638466/calculation	50000.0 ite.	ColumnOperator	[]
Run 3	TimeTrial: 0.41045224	29.47958252	6.64838434	36.70730908/calculation	50000.0 ite.	ColumnOperator	[]
Run 4	TimeTrial: 0.37967452	27.98394598	6.3926474	34.91495252/calculation	50000.0 ite.	ColumnOperator	[]
Run 5	TimeTrial: 0.35346734	27.69007386	6.20289318	34.39787556/calculation	50000.0 ite.	ColumnOperator	[]
Run 6	TimeTrial: 0.37319012	27.90951448	6.23438286	34.6706162/calculation	50000.0 ite.	ColumnOperator	[]
Run 7	TimeTrial: 0.3590824	26.99162072	6.18494232	33.68648468/calculation	50000.0 ite.	ColumnOperator	[]
Run 8	TimeTrial: 0.35533844	26.63953584	6.13542914	33.28009066/calculation	50000.0 ite.	ColumnOperator	[]
Run 9	TimeTrial: 0.35747358	26.8613118	6.20655488	33.58074044/calculation	50000.0 ite.	ColumnOperator	[]

Run 0	TimeTrial: 0.41600758	12.13486014	3.8182583	16.56913238/calculation	50000.0 ite.	BaseSingleOperator	[]
Run 1	TimeTrial: 0.21556538	11.62560862	3.72171942	15.73218808/calculation	50000.0 ite.	BaseSingleOperator	[]
Run 2	TimeTrial: 0.15663414	11.21613296	3.70233208	15.22494032/calculation	50000.0 ite.	BaseSingleOperator	[]
Run 3	TimeTrial: 0.17683894	12.65939916	3.83957374	16.83364216/calculation	50000.0 ite.	BaseSingleOperator	[]
Run 4	TimeTrial: 0.16401874	12.46954594	3.7763087	16.5608605/calculation	50000.0 ite.	BaseSingleOperator	[]
Run 5	TimeTrial: 0.17839254	12.22760226	3.81948232	16.38483204/calculation	50000.0 ite.	BaseSingleOperator	[]
Run 6	TimeTrial: 0.1647101	11.75324218	3.7241174	15.79404544/calculation	50000.0 ite.	BaseSingleOperator	[]
Run 7	TimeTrial: 0.16792772	11.00964618	3.79232506	15.12625314/calculation	50000.0 ite.	BaseSingleOperator	[]
Run 8	TimeTrial: 0.1926298	12.72196694	3.91817666	16.9991278/calculation	50000.0 ite.	BaseSingleOperator	[]
Run 9	TimeTrial: 0.2139492	12.9742752	4.05190914	17.40463786/calculation	50000.0 ite.	BaseSingleOperator	[]

TimeTrial: 433	4.33/calculation	FullMaster calculation no operator
TimeTrial: 194	1.94/calculation	Full calculation no operator
Run 0	TimeTrial: 0.4074847	38.6484565	6.93265446	46.22628088/calculation	50000.0 ite.	HaplotypeSwapSectionOperator	[5]
Run 1	TimeTrial: 0.22941214	34.69670638	6.374593	41.49219002/calculation	50000.0 ite.	HaplotypeSwapSectionOperator	[5]
Run 2	TimeTrial: 0.17901876	32.6449761	6.10610272	39.09047348/calculation	50000.0 ite.	HaplotypeSwapSectionOperator	[5]
Run 3	TimeTrial: 0.23603068	35.88257658	6.6243343	42.93206252/calculation	50000.0 ite.	HaplotypeSwapSectionOperator	[5]
Run 4	TimeTrial: 0.23365668	35.7456989	6.55165818	42.71487674/calculation	50000.0 ite.	HaplotypeSwapSectionOperator	[5]
Run 5	TimeTrial: 0.23418404	35.8877576	6.54972712	42.85458558/calculation	50000.0 ite.	HaplotypeSwapSectionOperator	[5]
Run 6	TimeTrial: 0.23364396	36.15824614	6.66463404	43.24176154/calculation	50000.0 ite.	HaplotypeSwapSectionOperator	[5]
Run 7	TimeTrial: 0.2474757	36.07406582	6.6105351	43.12473944/calculation	50000.0 ite.	HaplotypeSwapSectionOperator	[5]
Run 8	TimeTrial: 0.24822762	36.71654968	6.79185716	43.94988936/calculation	50000.0 ite.	HaplotypeSwapSectionOperator	[5]
Run 9	TimeTrial: 0.24830946	36.6983496	6.72501056	43.85674248/calculation	50000.0 ite.	HaplotypeSwapSectionOperator	[5]

Run 0	TimeTrial: 0.8593406	95.6122373	10.89572006	107.59898394/calculation	50000.0 ite.	BasesMultiOperator	[5]
Run 1	TimeTrial: 0.72309278	91.95091902	10.71005402	103.60769386/calculation	50000.0 ite.	BasesMultiOperator	[5]
Run 2	TimeTrial: 0.74440628	91.95233858	10.64581364	103.57091452/calculation	50000.0 ite.	BasesMultiOperator	[5]
Run 3	TimeTrial: 0.8185435	98.89591932	11.00317396	110.96230074/calculation	50000.0 ite.	BasesMultiOperator	[5]
Run 4	TimeTrial: 0.7843866	97.31593854	10.84113108	109.1776685/calculation	50000.0 ite.	BasesMultiOperator	[5]
Run 5	TimeTrial: 1.02011148	114.92207612	11.73383204	127.97256944/calculation	50000.0 ite.	BasesMultiOperator	[5]
Run 6	TimeTrial: 0.94567238	104.64419904	11.19189358	117.06038158/calculation	50000.0 ite.	BasesMultiOperator	[5]
Run 7	TimeTrial: 0.97568012	109.16105518	11.3715439	121.77248572/calculation	50000.0 ite.	BasesMultiOperator	[5]
Run 8	TimeTrial: 0.8310457	100.43465376	11.0025883	112.50818674/calculation	50000.0 ite.	BasesMultiOperator	[5]
Run 9	TimeTrial: 0.54303148	80.88046626	9.92586304	91.5475116/calculation	50000.0 ite.	BasesMultiOperator	[5]


*/