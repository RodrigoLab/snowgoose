package srp.likelihood.haplotypes;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.HashMap;

import javax.swing.text.TabableView;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.util.ArithmeticUtils;

import com.sun.org.glassfish.external.statistics.Stats;

import dr.math.BigDecimalUtils;
import dr.math.MathUtils;
import dr.math.Polynomial;
import dr.math.Polynomial.BigDouble;
import dr.math.distributions.NormalDistribution;
import srp.core.BigFunctions;
import srp.evolution.OperationType;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
//import srp.evolution.shortreads.AlignmentMapping;
import srp.evolution.shortreads.ShortReadMapping;
import srp.likelihood.AbstractShortReadsLikelihood;
import srp.likelihood.LikelihoodScaler;
//import srp.likelihood.spectrum.AbstractShortReadsSpectrumLikelihood;
import srp.likelihood.stateLikelihood.StateLikelihood;
//import java.util.BitSet;



public class ShortReadsHaplotypeLikelihood  extends AbstractShortReadsLikelihood {

	private static final long serialVersionUID = 7438385718398999755L;

	private static final boolean DEBUG = false;
	
//    public static final String SHORT_READ_LIKELIHOOD = "ShortReadHaplotypeLikelihood";
	public static final String NAME = "ShortReadHaplotypeLikelihood";

//	private final double MIN_LOG_LIKELIHOOD;

	

	
	protected HaplotypeModel alignmentModel;
	
	@Deprecated protected double[] allStateLogLikelihood;
	@Deprecated protected double[] allStoredStateLogLikelihood;

//	protected HaplotypeModel haplotypeModel;
	
	private HashMap<Integer, double[]> scaledLogBinomialDesnity;
	private HashMap<Integer, BigDecimal[]> bigDecBinomialDesnity;
	private int[][] allDists;
	private int[][] storedAllDists;
	
	@Deprecated protected StateLikelihood stateLikelihood;

//	private BigDecimal[] sumBigDecSrp;
//	private BigDecimal[] storedSumBigDecSrp;
//	private BigDecimal[] bigDecEachSrpLogLikelihood;
//	private BigDecimal[] storedBigDecEachSrpLogLikelihood;
//	private BigDecimal bigDecLikelihood;
//	private BigDecimal storedBigDecLikelihood;
	

	public ShortReadsHaplotypeLikelihood(HaplotypeModel haplotypeModel, ShortReadMapping srpMap){
		super(SHORT_READ_LIKELIHOOD, srpMap);
		this.alignmentModel = haplotypeModel;

		operationRecord = alignmentModel.getOperationRecord();
//		multiType = MultiType.Array;
		multiType = MultiType.BitSet;
		
//		type = MultiType.Hash;
//		type = MultiType.All;
//		distTypeCode = "flat";//"betaMean"  "betaMode" "gTest"
//		setDistType(distType);
//		MIN_LOG_LIKELIHOOD = 0;//stateLikelihood.caluclateStateLogLikelihood(SpectraParameter.MIN_FREQ);
		
		

		likelihoodKnown = false;
		
		addModel(this.alignmentModel);
		
		preprocessLikelihoodAlignmentMap();

//		calculateSrpLikelihoodFull();//TODO FIX this? shouldn't needed
		getLogLikelihood();
				
		storeEverything();
		
		for (int j = 0; j < sequenceCount; j++) {
			Haplotype haplotype = this.alignmentModel.getHaplotype(j);
			haplotype.storeState();
		}
		
	}
	
	

	private void preprocessLikelihoodAlignmentMap() {
//		makeDirty();


		
//		srpCount = srpMap.getSrpCount();
		sequenceCount = alignmentModel.getHaplotypeCount();
//		sequenceLength = alignmentModel.getHaplotypeLength();
//		operationRecord = alignmentModel.getOperationRecord();
//		this.srpSwitch = new boolean[srpCount];
//		this.allSrpPos = new HashSet<Integer>();
//		this.srpIndex = new int[srpCount];
		
//		logLikelihood = Double.NEGATIVE_INFINITY;
//		storedLogLikelihood = Double.NEGATIVE_INFINITY;

//		spectrumLogLikelihood = new double[srpCount*sequenceCount];
//		storedSpectrumLogLikelihood = new double[srpCount*sequenceCount];
//
//		scaledSpectrumLogLikelihood = new double[srpCount*sequenceCount];
//		storedScaledSpectrumLogLikelihood = new double[srpCount*sequenceCount];
//		
//		sumScaledSrpLogLikelihood = new double[srpCount];
//		storedSumSrpLogLikelihood = new double[srpCount];
//
//		
//		eachSrpLogLikelihood = new double[srpCount];
//		storedEachSrpLogLikelihood = new double[srpCount];

//		allStateLogLikelihood = new double[sequenceLength*STATE_COUNT];
//		allStoredStateLogLikelihood = new double[sequenceLength*STATE_COUNT];
		
//		mapToSrpArray = srpMap.getMapToSrpArray();
//		
//		
//		String[] srpArray = srpMap.getSrpArray();
//		
//		allSrpState2D = new int[srpArray.length][sequenceLength];
//		allSrpChar2D = new char[srpArray.length][sequenceLength];
//		
//		for (int i = 0; i < srpArray.length; i++) {
//			String srp = srpArray[i];
//			for (int j = 0; j < sequenceLength; j++) {
//				allSrpState2D[i][j] = getStateAtK(srp, j);
//				allSrpChar2D[i][j] = srp.charAt(j);
//			}
//		}
//		
//		
//		bitSet = new BitSet(srpCount);
//		sumBigDecSrp = new BigDecimal[srpCount];
//		storedSumBigDecSrp = new BigDecimal[srpCount];
//		bigDecEachSrpLogLikelihood = new BigDecimal[srpCount];
//		storedBigDecEachSrpLogLikelihood = new BigDecimal[srpCount];
//		bigDecLikelihood = new BigDecimal("0");
//		storedBigDecLikelihood = new BigDecimal("0");
		
		
		scaledLogBinomialDesnity = new HashMap<Integer, double[]>();
		bigDecBinomialDesnity = new HashMap<Integer, BigDecimal[]>();
//		scaledLogBinomialDesnity = new double[]
		allDists = new int[srpCount][sequenceCount];
		storedAllDists = new int[srpCount][sequenceCount];

		double[] logFractorial = new double[500];//Change 500 later TODO:
		for (int i = 0; i < logFractorial.length; i++) {
			logFractorial[i] = ArithmeticUtils.factorialLog(i);
		}
		
		double min = 100;
		BigDecimal BigE = new BigDecimal("0.107");
		BigDecimal BigN = new BigDecimal("0.893");
		for (int s = 0; s < srpCount; s++) {
			int srLength = srpMap.getSrpLength(s);//ength();//srp.length();
			int srLengthPlusOne = srLength+1;
//			maxDist = Math.max(maxDist, srLength1);
			double[] binomD = new double[srLengthPlusOne];
			BigDecimal[] bigBinom = new BigDecimal[srLengthPlusOne];
			double[] logBinomD = new double[srLengthPlusOne];
			double[] scaledBinomD = new double[srLengthPlusOne];
			double sum=0;
			
			
			double lambda = srLength * ERROR_RATE;
			double logLambda = Math.log(lambda);
			for (int i = 0; i < logBinomD.length; i++) {
//				BigDecimal t1 = new BigDecimal(String.valueOf(i));
//				BigDecimal t2 = new BigDecimal(String.valueOf(LOG_ERROR_RATE));
//				BigDecimal t3 = new BigDecimal(String.valueOf(srLength-i));
//				BigDecimal t4 = new BigDecimal(String.valueOf(LOG_ONE_MINUS_ERROR_RATE));
//				BigDecimal t5 = t1.multiply(t2).add(t3.multiply(t4));
////				*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE)
				logBinomD[i] = i*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE;
				
//				logBinomD[i] = i* logLambda - logFractorial[i] - lambda;
//				System.out.println(i +"\t"+  srLength + "\t"+ lambda +"\t"+ logBinomD[i]);
//				logBinomD[i] = (i*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE) /10;
				logBinomD[i] = ArithmeticUtils.binomialCoefficientLog(srLength, i)+i*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE;
//				System.out.println( (binomD[i]==Math.exp(logBinomD[i])) +"\t"+ binomD[i] +"\t"+ Math.exp(logBinomD[i]) +"\t"+ Math.log(binomD[i]) +"\t"+ logBinomD[i] );
				bigBinom[i] = new BigDecimal("0");
//				bigBinom[i].setScale(1000);
				BigDecimal temp = BigE.pow(i).multiply(BigN.pow(srLength-i));
//				System.out.println(bigBinom[i].scale() +"\t"+ bigBinom[i].toPlainString());
				bigBinom[i] = bigBinom[i].add(temp);
//				BigDecimal temp2 = temp.add(bigBinom[i]);
//				System.out.println(temp2.doubleValue() +"\t"+ bigBinom[i].toPlainString());
//				bigBinom[i] = new BigDecimal(multiply.toPlainString()); // better than bigdb2 and bigdb3
//				BigDecimal bigbd2 = BigFunctions.exp(new BigDecimal(Double.toString(logBinomD[i])), 500);
//				BigDecimal bigbd3 = BigFunctions.exp(t5, 500);
				sum += Math.exp(logBinomD[i]);
//				binomD[i] = Math.exp(logBinomD[i]);
//				if(min > binomD[i]){
//					min = binomD[i];
//				}
//				bigBinom[i] = bigbd;//BigDecimal.valueOf(binomD[i]);
//				bigBinom[i].setScale(50);
//				System.out.println(srLength +"\t"+ i);
//				System.out.println(binomD[i] +"\n"+ bigBinom[i].toPlainString());
//				bigbd = bigbd.setScale(50, BigDecimal.ROUND_HALF_EVEN);
//				System.out.println(bigBinom[i].scale() +"\t"+  bigBinom[i].toPlainString());
//				System.out.println(bigbd.scale() +"\t"+  bigbd.toPlainString());
//				System.out.println(bigbd2.scale() +"\t"+  bigbd2.toPlainString());
//				System.out.println(bigbd3.scale() +"\t"+  bigbd3.toPlainString());
//				System.out.println();
//				bigbd.s
//				BinomialDistribution
//				NormalDistribution
//				logBinomD[i] = i*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE;
				
//				System.out.println(logBinomD[i] +"\t"+ scaledBinomD[i]);
//				scaledBinomD[i] = i/100000.0;
//				scaledBinomD[i] = Math.pow(LOG_ERROR_RATE, i)*Math.pow(LOG_ONE_MINUS_ERROR_RATE, (srLength-i));
				
				scaledBinomD[i] = liS.scale(logBinomD[i]);
			}
			double logSum = Math.log(sum);
//			double sum2 = 0;
			for (int i = 0; i < logBinomD.length; i++) {
				logBinomD[i] = logBinomD[i] -logSum;
//				sum2 += Math.exp(logBinomD[i]);
				scaledBinomD[i] = liS.scale(logBinomD[i]);
			}
//			System.out.println(srLength +"\t"+ sum2 +"\t"+ Arrays.toString(logBinomD)) ;
	//		System.out.println(Arrays.toString(logBinomD));
//			logBinomialDesnity.put(srLength, logBinomD);
			scaledLogBinomialDesnity.put(srLength, scaledBinomD);
			bigDecBinomialDesnity.put(srLength, bigBinom);
		}
//		System.out.println(min);
	}


    
	protected double calculateSrpLikelihoodFull() {

//		System.out.println("calculateSrpLikelihood_Full");
		double logLikelihood = 0;
//		bigDecLikelihood = new BigDecimal("0");
		for (int i = 0; i < srpCount; i++) {
			
			String srp = srpMap.getSrpFragment(i);//srpArray[i]
			int start = srpMap.getSrpStart(i);
			int end = srpMap.getSrpEnd(i);
			
			double[] logPD = scaledLogBinomialDesnity.get(srpMap.getSrpLength(i));
			liS.reset();
			sumScaledSrpLogLikelihood[i] = 0;
			for (int j = 0; j < sequenceCount; j++) {

				int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(j));
				allDists[i][j]=dist;
				liS.add(logPD[dist]);
//				System.out.println(sumBigDecSrp[i]);
//				sumScaledSrpLogLikelihood[i] += logPD[allDists[i][j]];
//				storedAllLog2D[i][j] = logPD[dist];
//				liS.scaleLogProb(logPD[dist]);
			}	
			sumScaledSrpLogLikelihood[i] = liS.getSumScaledLikelihood();
			eachSrpLogLikelihood[i] = liS.getLogLikelihood() ;
			
			eachSrpLogLikelihood[i] = liS.getLogLikelihood();
//			bigDecEachSrpLogLikelihood[i] = BigFunctions.ln(sumBigDecSrp[i], 300);
//			eachSrpLogLikelihood[i] = bigDecEachSrpLogLikelihood[i].doubleValue();
//			System.out.println(eachSrpLogLikelihood[i] +"\t"+ sumBigDecSrp[i].toPlainString());
//			System.out.println("inFull: "+sumScaledSrpLogLikelihood[i] +"\t"+ eachSrpLogLikelihood[i]);
//			bigDecLikelihood = bigDecLikelihood.add(bigDecEachSrpLogLikelihood[i]);
			logLikelihood += eachSrpLogLikelihood[i];
		}
//		System.out.println(logLikelihood +"\t"+ bigDecLikelihood.toPlainString());
		return logLikelihood;
	}

	static String cheat[] = new String[]{
		"AGTTAAGAGGCACAACTTTGGTGGATGTCAGTTGAGGTTGCTACTACACAAAAAGTACATACCTTACGGAACGTTCTCAATCACTGTGAGACGTTTGCCATGGTACACTCTGCGTCCACT",
		"AGTTAAGAGGCACAACTTTGGTGGATGTCAGTTGAGGTTGCTACTACACAAAAAGTACATACCTTACGGCACGTTCTCAATCACTGTGAGACGTTTGCCATGGTACACTCTGCGTCCACT",
		"AGTTAAGAGGCACAACTTTGGTGGATGTCAGTTGAGGTTGCTACTACACAAAAAGTACATACCTTACGGAACGTTCTCAATCACTGTGAGACGTTTGCCATGGTACACTCTGCGTCCACT",
		"AGCTAAGAGACACAACTATGGTGGATGTCAGTCAAGCCTGCAACTACATAGAAAGTGAATAAGTTACTGAATGCTCTCATTCACTGAGAGACGTTTACCATAATATACTAGGCGTCGACG",
		"AGCTAAGAGACACAACTATGGTGGATGTCAGTCAAGCCTGCAACTACATAGAAAGTGAATAAGCTACTGAATGCTCTCATTCACTGAGAGACGTTTACCATAATATACTAGGCGTCGACG",
		"AGCTAAGAGACACAACTATGGTGGATGTCAGTCAAGCCTGCAACTACATAGAAAGTGAATAAGTTACTGAATGCTCTCATTCACTGAGAGACGTTTACCATAATATACTTGGCGTCGACG",
		"AGCTAAGAGACACAACTATGGTGGATGTCAGTCAAGCCTCCTACTACATAGAAAGTGAATAAGTTACTGAATGCTCTCATTCACTGAGAGACGTTTGCCTTAATAGACTTGGCGTCGACG",
	};
	
	private int calculateDeltaDist(char srpChar, char newChar, char oldChar) {
			int deltaDist = 0;
			
			if (srpChar==newChar){
				deltaDist = -1;
			}
			else if(srpChar==oldChar){
				deltaDist = 1;
			}
	//		System.out.println("\tCalDelta: "+srpChar +"\t"+ newChar +"\t"+ oldChar +"\t"+ deltaDist);
			return deltaDist;
		}



	protected double calculateSrpLikelihoodSingle() {

		int j = operationRecord.getSpectrumIndex(); 
		int k = operationRecord.getSingleIndex();//AllSiteIndexs()[0];
		
		double currentLogLikelihood = getStoredLogLikelihood();

		Haplotype haplotype = alignmentModel.getHaplotype(j);
		char oldChar = haplotype.getStoredChar(k);
		char newChar = haplotype.getChar(k);
		
//		System.out.print("\t"+ currentLogLikelihood +"\t" );
		if(newChar!= oldChar){ 
			for (int i : mapToSrpArray[k]){
				char srpChar = allSrpChar2D[i][k];
				int deltaDist = calculateDeltaDist(srpChar, newChar, oldChar);
				
				if (deltaDist != 0) {
					int srpLength = srpMap.getSrpLength(i);
//					BigDecimal[] tempBD = bigDecBinomialDesnity.get(srpLength);
					double[] logPD = scaledLogBinomialDesnity.get(srpLength);
					
					double sum = 0;
					
//					double sum = 0;
//					BigDecimal bd = new BigDecimal("0");
//					bd.setScale(100);
					for (int l = 0; l < sequenceCount; l++) {
						sum += logPD[allDists[i][l]];
//						bd = bd.add(BigDecimal.valueOf( logPD[allDists[i][l]] ));
					}
					if(deltaDist == -1){
//						System.out.println(deltaDist +"\t"+  srpChar +"\t"+ newChar +"\t"+ oldChar +"\t"+ cheat[j].charAt(k));
						System.out.println("Old\t"+currentLogLikelihood +"\t"+ eachSrpLogLikelihood[i] );
					}
//					sumBigDecSrp[i] = sumBigDecSrp[i].subtract(tempBD[allDists[i][j]]);
					
					allDists[i][j] += deltaDist;
					
//					sumBigDecSrp[i] = sumBigDecSrp[i].add(tempBD[allDists[i][j]]);
					sum = 0;
					for (int l = 0; l < sequenceCount; l++) {
						sum += logPD[allDists[i][l]];
//						bd = bd.add(BigDecimal.valueOf( logPD[allDists[i][l]] ));
					}
//					bigDecLikelihood = bigDecLikelihood.subtract(bigDecEachSrpLogLikelihood[i]);
					currentLogLikelihood -= eachSrpLogLikelihood[i];
//					BigDecimal bd = BigFunctions.ln(sumBigDecSrp[i], sumBigDecSrp[i].scale());
					eachSrpLogLikelihood[i] = LikelihoodScaler.getLogLikelihood(sum, LOG_C);
//					eachSrpLogLikelihood[i]
//					bigDecEachSrpLogLikelihood[i] = BigFunctions.ln(sumBigDecSrp[i], 300);
//					eachSrpLogLikelihood[i] = bigDecEachSrpLogLikelihood[i].doubleValue();
//					bigDecLikelihood = bigDecLikelihood.add(bigDecEachSrpLogLikelihood[i]);
					
					
//					eachSrpLogLikelihood[i] = bd.doubleValue();
//					eachSrpLogLikelihood[i] = liS.getLogLikelihood();//LiS
					currentLogLikelihood += eachSrpLogLikelihood[i];
					
//					currentLogLikelihood = updateEachSrpAtI(i, currentLogLikelihood, (allDists[i][j]-deltaDist), (allDists[i][j]));
					if(deltaDist == -1){
						double newSum = 0;
						double allLogSum = 0;
//						BigDecimal newBd = new BigDecimal("0");
//						newBd.setScale(100);
//						int srpLength = srpMap.getSrpLength(i);
//						double[] logPD = scaledLogBinomialDesnity.get(srpLength);
						for (int l = 0; l < sequenceCount; l++) {
//							newSum += allDists[i][l];
							newSum += logPD[allDists[i][l]];
//							BigDecimal tempB = BigDecimal.valueOf( logPD[allDists[i][l]] );
//							
//							newBd = newBd.add(tempB);
//							
						}
						System.out.println("new\t"+currentLogLikelihood +"\t"+ eachSrpLogLikelihood[i] ); 
//						sumBigDecSrp[i].doubleValue() +"\t"+  sumBigDecSrp[i].toPlainString());
//						System.out.println(LikelihoodScaler.getLogLikelihood(sum, LOG_C) +"\t"+ bigDecEachSrpLogLikelihood[i].doubleValue() +"\t"+ bigDecEachSrpLogLikelihood[i].toPlainString());
//						System.out.println(currentLogLikelihood +"\t"+ bigDecLikelihood.doubleValue() +"\t"+ bigDecLikelihood.toPlainString());
//						System.out.println("===\t"+ newBd.doubleValue());
						double log2 = newSum*LOG_ERROR_RATE+(srpLength*sequenceCount-newSum)*LOG_ONE_MINUS_ERROR_RATE;
//						System.out.println(sum +"\t"+ bd.toString() +"\n"+ newSum +"\t"+ newBd.toString() );  
//						System.out.println(liS.getLogLikelihood(bd.doubleValue(), LOG_C) +"\t"+ liS.getLogLikelihood(newBd.doubleValue(), LOG_C) );
//						System.out.println(liS.getLogLikelihood(Double.valueOf(bd.toPlainString()), LOG_C) +"\t"+ liS.getLogLikelihood(Double.valueOf(newBd.toPlainString()), LOG_C) );
//						System.out.println(BigFunctions.ln(bd, 100).toPlainString() +"\n"+ BigFunctions.ln(newBd, 100).toPlainString());
//						System.out.println((BigFunctions.ln(bd, 100).doubleValue())-LOG_C +"\n"+ (BigFunctions.ln(newBd, 100).doubleValue()-LOG_C));
						System.out.println(allDists[i][j] +"\t"+  logPD[allDists[i][j]] +"\t"+ logPD[ (allDists[i][j]-deltaDist  )]);
						System.out.println(sumScaledSrpLogLikelihood[i] +"\t"+  storedSumScaledSrpLogLikelihood[i] +"\t");//storedBigDecLikelihood.toPlainString());
						
						System.out.println();
						System.out.println(LOG_C);
					}
				}

			}
		}
//		currentLogLikelihood = bigDecLikelihood.doubleValue();
//		System.out.println(currentLogLikelihood);
		return currentLogLikelihood;
	}

	protected double calculateSrpLikelihoodSingle3() {

		int j = operationRecord.getSpectrumIndex(); 
		int k = operationRecord.getSingleIndex();//AllSiteIndexs()[0];
		
		double currentLogLikelihood = getStoredLogLikelihood();

		Haplotype haplotype = alignmentModel.getHaplotype(j);
		char oldChar = haplotype.getStoredChar(k);
		char newChar = haplotype.getChar(k);
		
//		System.out.print("\t"+ currentLogLikelihood +"\t" );
		if(newChar!= oldChar){ 
			for (int i : mapToSrpArray[k]){
				char srpChar = allSrpChar2D[i][k];
				int deltaDist = calculateDeltaDist(srpChar, newChar, oldChar);
				
				if (deltaDist != 0) {
					int srpLength = srpMap.getSrpLength(i);
					double[] logPD = scaledLogBinomialDesnity.get(srpLength);
					
					double sum = 0;
					BigDecimal bd = new BigDecimal("0");
					bd.setScale(100);
					for (int l = 0; l < sequenceCount; l++) {
						sum += logPD[allDists[i][l]];
						bd = bd.add(BigDecimal.valueOf( logPD[allDists[i][l]] ));
					}
					
					allDists[i][j] += deltaDist;
					if(deltaDist == -1){
						System.out.println(deltaDist +"\t"+  srpChar +"\t"+ newChar +"\t"+ oldChar +"\t"+ cheat[j].charAt(k));
						System.out.println(currentLogLikelihood +"\t"+ eachSrpLogLikelihood[i]);
					}
					currentLogLikelihood = updateEachSrpAtI(i, currentLogLikelihood, (allDists[i][j]-deltaDist), (allDists[i][j]));
					if(deltaDist == -1){
						double newSum = 0;
						double allLogSum = 0;
						BigDecimal newBd = new BigDecimal("0");
						newBd.setScale(100);
//						int srpLength = srpMap.getSrpLength(i);
//						double[] logPD = scaledLogBinomialDesnity.get(srpLength);
						for (int l = 0; l < sequenceCount; l++) {
//							newSum += allDists[i][l];
							newSum += logPD[allDists[i][l]];
							BigDecimal tempB = BigDecimal.valueOf( logPD[allDists[i][l]] );
							
							newBd = newBd.add(tempB);
							
						}
//						System.out.println("===\t"+ newBd.doubleValue());
//						double log2 = newSum*LOG_ERROR_RATE+(srpLength*sequenceCount-newSum)*LOG_ONE_MINUS_ERROR_RATE;
						System.out.println(sum +"\t"+ bd.toString() +"\n"+ newSum +"\t"+ newBd.toString() );  
						System.out.println(liS.getLogLikelihood(bd.doubleValue(), LOG_C) +"\t"+ liS.getLogLikelihood(newBd.doubleValue(), LOG_C) );
						System.out.println(liS.getLogLikelihood(Double.valueOf(bd.toPlainString()), LOG_C) +"\t"+ liS.getLogLikelihood(Double.valueOf(newBd.toPlainString()), LOG_C) );
						System.out.println(BigFunctions.ln(bd, 100).toPlainString() +"\n"+ BigFunctions.ln(newBd, 100).toPlainString());
						System.out.println((BigFunctions.ln(bd, 100).doubleValue())-LOG_C +"\n"+ (BigFunctions.ln(newBd, 100).doubleValue()-LOG_C));
						System.out.println(allDists[i][j] +"\t"+  logPD[allDists[i][j]] +"\t"+ logPD[ (allDists[i][j]-deltaDist  )]);
						System.out.println(sumScaledSrpLogLikelihood[i] +"\t"+  storedSumScaledSrpLogLikelihood[i]);
						System.out.println(currentLogLikelihood +"\t"+ eachSrpLogLikelihood[i]);
						System.out.println();
//						System.out.println(LOG_C);
					}
				}

			}
		}
//		System.out.println(currentLogLikelihood);
		return currentLogLikelihood;
	}
	

	
	protected double calculateSrpLikelihoodSingle2() {


//		OperationRecord record = alignmentModel.getOperationRecord();
		int j = operationRecord.getSpectrumIndex(); 
		int k = operationRecord.getSingleIndex();//AllSiteIndexs()[0];
		double currentLogLikelihood = getStoredLogLikelihood();
		
//		System.out.println("StartSingle" +"\t"+ j +"\t"+ k +"\t"+ currentLogLikelihood);
//		SpectraParameter spectra = haplotypeModel.getHaplotype(j).getSpectra(k);
//		ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);

//		stateLikelihood.calculateStatesLogLikelihood2D(spectra, allStateLogLikelihood2D[0]);
//		stateLikelihood.calculateStoredStatesLogLikelihood2D(spectra, allStoredStateLogLikelihood2D[0]);

//		stateLikelihood.calculateStatesLogLikelihood(spectra, 0, allStateLogLikelihood);
//		stateLikelihood.calculateStoredStatesLogLikelihood(spectra, 0, allStoredStateLogLikelihood);

//		for (int i : mapToSrp) {
		Haplotype haplotype = alignmentModel.getHaplotype(j);
		char oldChar = haplotype.getStoredChar(k);
		char newChar = haplotype.getChar(k);

		
		if(newChar!= oldChar){ // if(newChar!= oldChar && isHapEqualNew)
			for (int i : mapToSrpArray[k]){
				int srpLength = srpMap.getSrpLength(i);
	//			int state = getStateAtK(fullSrp, k);
	//			int state = allSrpState2D[i][k];
				
				char srpChar = allSrpChar2D[i][k];
				int deltaDist = 0;
				
				if (srpChar==newChar){
					deltaDist = -1;
				}
				else if(srpChar==oldChar){
					deltaDist = 1;
				}
				
//				
//				String srp = srpMap.getSrpFragment(i);//srpArray[i]
//				int start = srpMap.getSrpStart(i);
//				int end = srpMap.getSrpEnd(i);
//				double[] logPD = scaledLogBinomialDesnity.get(srpMap.getSrpLength(i));
//				liS.reset();
//				double addLogDist = 0;
//				for (int hj = 0; hj < sequenceCount; hj++) {
//					
//					int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(hj).getSequenceString());
////					System.out.println(hj +" "+ j +":\t"+ dist +"\t"+ storedAllDists[i][hj] +"\t"+ (dist == storedAllDists[i][hj]));
//					liS.add(logPD[dist]);
//					addLogDist += logPD[dist];
//					if( storedAllLog2D[i][hj] != logPD[dist]){
//						System.out.println("ChangeHap: "+(j==hj) +"\t"+ i +"\t"+ j +"\t"+ hj +"\tStored: "+storedAllLog2D[i][hj] +"\t"+  logPD[dist] +"\t\t"+ srpMap.getSrpLength(i));
//					}
//					if(dist != storedAllDists[i][hj]){
//						System.out.println("ChangeHap: "+(j==hj) +"\t"+ i +"\t"+ j +"\t"+ hj +"\tNewDist:"+ dist +"\t"+ storedAllDists[i][hj] );
//					}
//				}
//				double e1 = liS.getSumScaledLikelihood(); 
//				double e2 = liS.getLogLikelihood();
				
				if (deltaDist != 0) {
					
//					int srpChar = srp.getFullSrpCharAt(swapPos);
//					int deltaDist = calculateDeltaDist(srpChar, newChar, oldChar);

					double[] logPD = scaledLogBinomialDesnity.get(srpLength);
					int oldDist = storedAllDists[i][j];
					int newDist = storedAllDists[i][j] + deltaDist;
					allDists[i][j] = newDist;
//					storedAllLog2D[i][j] = logPD[newDist];
//					
//					double oldPD = logPD[oldDist];
//					double newPD = logPD[newDist];
//					
//					liS.reset();		
//					for (int hj = 0; hj < sequenceCount ; hj++) {
//						if(j == hj){
//							liS.add(logPD[storedAllDists[i][hj]]);
//						}
//						else{
//							liS.add(logPD[allDists[i][hj]]);
//						}
//					}
//					double sumSrp = liS.getSumScaledLikelihood();
//					double eachSrp = liS.getLogLikelihood();
//					if(eachSrp != eachSrpLogLikelihood[i]){
//						System.out.println("check stored??" +"\t"+ eachSrp +"\t"+ eachSrpLogLikelihood[i]);
//					}
//					if(sumSrp != sumScaledSrpLogLikelihood[i]){
//						System.out.println("check stored sum srp??" +"\t"+ sumSrp +"\t"+ sumScaledSrpLogLikelihood[i]);
//					}
					
					currentLogLikelihood -= eachSrpLogLikelihood[i];
//					sumScaledSrpLogLikelihood[i] -= logPD[oldDist];//scaledSpectrumLogLikelihood[offset];

//						spectrumLogLikelihood[offset] -= storedStateLn; 
//						spectrumLogLikelihood[offset] += stateLn;
//						scaledSpectrumLogLikelihood[offset] = LikelihoodScaler.scale(spectrumLogLikelihood[offset], LOG_C);
					
					sumScaledSrpLogLikelihood[i] = storedSumScaledSrpLogLikelihood[i] - (logPD[oldDist] - logPD[newDist]);//scaledSpectrumLogLikelihood[offset];
					double srpSum = 0;
					double neg = sumScaledSrpLogLikelihood[i];
					for (int hj = 0; hj < sequenceCount; hj++) {
						srpSum += logPD[allDists[i][hj]];
//						if(hj == j){
//							neg -= logPD[storedAllDists[i][hj]];
//						}
//						else{
							neg -= logPD[allDists[i][hj]];
//						}
					}
					 srpSum = 0;
					double oldSum = 0;
					for (int hj = 0; hj < sequenceCount; hj++) {
						srpSum += logPD[allDists[i][hj]];
						oldSum += logPD[storedAllDists[i][hj]];
						if(allDists[i][hj]!=storedAllDists[i][hj]){
//							if(j != hj){
								System.out.println(i +"\t"+  j +"\t"+ hj +"\t"+ allDists[i][hj] +"\t"+ storedAllDists[i][hj] );
//								System.exit(-1);
//							}
						}
					}
					double x1 = (srpSum - logPD[allDists[i][j]]); 
					double x2 = (oldSum - logPD[storedAllDists[i][j]]);
					if(x1!=x2){
						System.out.println(srpSum +"\t"+ oldSum +"\t"+ (oldSum ==storedSumScaledSrpLogLikelihood[i]));
						System.out.println(allDists[i][j] +"\t"+ storedAllDists[i][j]);
						System.out.println( x1 +"\t"+ x2); 
						
						System.exit(-1);
					}
					String srp = srpMap.getSrpFragment(i);//srpArray[i]
					int start = srpMap.getSrpStart(i);
					int end = srpMap.getSrpEnd(i);
					int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(j));
					int dist2 = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(j).getStoredSequenceString());
					if(allDists[i][j] != dist){
						System.out.println("Start: "+start +"\t"+ end );
						System.out.println(srp);
						System.out.println(alignmentModel.getHaplotype(j).getSequenceString().substring(start,end));
						System.out.println(alignmentModel.getHaplotype(j).getStoredSequenceString().substring(start, end));
//						for (int k : siteIndexs) {
							 oldChar = haplotype.getStoredChar(k);
							 newChar = haplotype.getChar(k);
//							 srpChar = srpChars[k];
							 deltaDist = 0;
							if(newChar!= oldChar){ 
								if (srpChar==newChar){
									deltaDist = -1;
								}
								else if(srpChar==oldChar){
									deltaDist = 1;
								}
							}
							System.out.println(k +"\t"+ srpChar +"\t"+ newChar +"\t"+ oldChar +"\t"+ deltaDist);
//						}
						System.out.println(i +"\t"+ j +"\t"+ dist +"\t"+ dist2 +"\t"+ allDists[i][j]+"\t"+storedAllDists[i][j] +"\t");
						System.exit(-1);
					}
//					for (int hj = 0; hj < sequenceCount; hj++) {
//						int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(hj).getSequenceString());
//						if(allDists[i][hj] != dist){
//							System.out.println(dist +"\t"+ allDists[i][hj]+"\t"+ i +"\t"+ hj);
//						}
//					}
					double oldUpdate = oldSum -  - logPD[oldDist] + logPD[newDist];
					sumScaledSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i] - logPD[oldDist] + logPD[newDist];//scaledSpectrumLogLikelihood[offset];
					if(srpSum != sumScaledSrpLogLikelihood[i]){
						double d1 = srpSum-logPD[newDist];
						double d2 = storedSumScaledSrpLogLikelihood[i] - logPD[oldDist];
						System.out.println((d1==d2)+"\t"+ d1 +"\t"+ d2 +"\t");
						System.out.println((oldSum==storedSumScaledSrpLogLikelihood[i]) +"\t"+ (oldUpdate==srpSum) +"\t"+ (oldUpdate == sumScaledSrpLogLikelihood[i]) );
						System.out.println(srpSum +"\t"+ (srpSum-sumScaledSrpLogLikelihood[i]));
					}
//					if(srpSum != sumScaledSrpLogLikelihood[i]){
//						System.out.println("different sum: "+"\t"+ srpSum +"\t"+ sumScaledSrpLogLikelihood[i] );
//						System.out.println("\t\t"+ (neg==0) +"\t"+ neg +"\t"+ (logPD[oldDist] - logPD[newDist]));
//						System.exit(-1);
//					}
//					eachSrpLogLikelihood[i] = LikelihoodScaler.getLogLikelihood(sumScaledSrpLogLikelihood[i], LOG_C);
					eachSrpLogLikelihood[i] = LikelihoodScaler.getLogLikelihood(srpSum, LOG_C);
					currentLogLikelihood += eachSrpLogLikelihood[i];
					
//					sumSrp = sumSrp - oldPD + newPD;
//					liS.reset();
//					liS.add(sumSrp);
//					eachSrp = liS.getLogLikelihood();
//					
//					if(eachSrp != eachSrpLogLikelihood[i]){
//						System.out.println("storing srp error??" +"\t"+ eachSrp +"\t"+ eachSrpLogLikelihood[i]);
//					}
//					
//					if(sumSrp != sumScaledSrpLogLikelihood[i]){
//						System.out.println("storing srp error??" +"\t"+ sumSrp +"\t"+ sumScaledSrpLogLikelihood[i]);
//					}
//					System.out.println("\tdeltaDist: " + deltaDist +"\tnewDist: "+ newDist +"\toldDist: "+ oldDist);
//				
//					
//					
//					System.out.println((e1==addLogDist) +" liS working?" );
//					
//					double d1 = addLogDist - newPD;
//					double d2 = storedSumScaledSrpLogLikelihood[i] - oldPD;
//					double d3 = sumScaledSrpLogLikelihood[i] - newPD;
//					System.out.println((d2==d3) +" storing working?" +"\t"+ d2 +"\t"+ d3);
//					
//					System.out.println((e1 == sumScaledSrpLogLikelihood[i]) +"\t"+ e1 +"\t"+ sumScaledSrpLogLikelihood[i]);
////					System.out.println("\t"+ newPD +"\t"+ oldPD);
//					System.out.println((d1==d2) +"\t"+ d1 +"\t"+ d2);
//					
//					
//					if(e1 != 	sumScaledSrpLogLikelihood[i]){
//						System.out.println(i +"\t"+ eachSrpLogLikelihood[i] + "\t" + e2 +"\t"+ 
//								(eachSrpLogLikelihood[i] -e2) +"\t"+ 
//								e1 +"\t"+  
//								sumScaledSrpLogLikelihood[i]
//								);
//						System.exit(-1);
//					}

				}
				
			}
		}
//		//REMOVE
//		double trueL = calculateSrpLikelihoodFullMaster();
//		if((trueL - currentLogLikelihood)> 1e-8) {
//			System.out.println("DEBUG "+trueL +"\t"+ currentLogLikelihood);
//			
//			for (int i: mapToSrpArray[k]) {
//				
//				String srp = srpMap.getSrpFragment(i);//srpArray[i]
//				int start = srpMap.getSrpStart(i);
//				int end = srpMap.getSrpEnd(i);
//				
//				double[] logPD = scaledLogBinomialDesnity.get(srpMap.getSrpLength(i));
//
//				liS.reset();
//				double oldPD = 0;
//				double newPD = -1;
//				for (int hj = 0; hj < sequenceCount; hj++) {
//
//					int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(hj).getSequenceString());
////					System.out.println(hj +" "+ j +":\t"+ dist +"\t"+ storedAllDists[i][hj] +"\t"+ (dist == storedAllDists[i][hj]));
//					liS.add(logPD[dist]);
//					
//					if(hj == j){
//						int srpLength = srpMap.getSrpLength(i);
//						char srpChar = allSrpChar2D[i][k];
//						int deltaDist = 0;
//						if (srpChar==newChar){
//							deltaDist = -1;
//						}
//						else if(srpChar==oldChar){
//							deltaDist = 1;
//						}
//						if (deltaDist != 0) {
//		
//							int srpIndex = i; int hapIndex = j; 
//			
//							logPD = scaledLogBinomialDesnity.get(srpLength);
//							int oldDist = storedAllDists[srpIndex][hapIndex];
//							int newDist = storedAllDists[srpIndex][hapIndex] + deltaDist;
////							System.out.println();
//							oldPD = logPD[oldDist];
//							newPD = logPD[newDist];
//							System.out.println(hj + " " + j + ":\t" + dist +"\t"+ newDist
//									+ "\t" + storedAllDists[i][hj] + "\t" +
//									logPD[oldDist] +"\t"+ logPD[newDist] +"\t"+ 
//									(dist == storedAllDists[i][hj]));
//						}
//					}
////					liS.scaleLogProb(logPD[dist]);
//				}	
//							
//			}
//		}
		return currentLogLikelihood;
	}
	
	
	protected double calculateSrpLikelihoodMulti() {
		
		int[] siteIndexs = operationRecord.getAllSiteIndexs();
		int j= operationRecord.getSpectrumIndex(); 
		Haplotype haplotype = alignmentModel.getHaplotype(j);
		
		int multihere;
		
		double currentLogLikelihood = getStoredLogLikelihood();
		if(multiType == MultiType.BitSet){
			recalculateBitSet(siteIndexs);
			
			int count = 0;
			for (int i = bitSet.nextSetBit(0); i >= 0; i = bitSet.nextSetBit(i+1)) {
				srpIndex[count++] = i;
//				System.out.print(i +"\t"+ currentLogLikelihood);
				currentLogLikelihood = updateLikelihoodAtIJ(i, j,
						siteIndexs, haplotype, currentLogLikelihood);

			}
//			System.out.print(currentLogLikelihood);
//			double masterLog = calculateSrpLikelihoodFullMaster();
//			System.out.println("\t"+masterLog);
			srpIndexCount = count;

			
		}
		else if(multiType==MultiType.Array){
			
			recalculateArray(siteIndexs);

//			for (int i = 0; i < srpSwitch.length; i++) {
//				if(srpSwitch[i]){
//
//					
//				}
//			}
			throw new IllegalArgumentException("Not yet implemented");
		}
		else if(multiType==MultiType.Hash){
			recalculateHashSet(siteIndexs);

//			for (int i : allSrpPos) {
//				
//			}
			throw new IllegalArgumentException("Not yet implemented");
		}
		
		
		

		return currentLogLikelihood;

	}


	protected double calculateSrpLikelihoodColumn() {

//		OperationRecord record = alignmentModel.getOperationRecord();
		int k = operationRecord.getSingleIndex();

		int[] allSpectrumIndexs = operationRecord.getAllSpectrumIndexs();
		double currentLogLikelihood = getStoredLogLikelihood();
		
//		for (int j : allSpectrumIndexs) {
//			SpectraParameter spectra = haplotypeModel.getHaplotype(j).getSpectra(k);
//			
//			int kOffset = j*STATE_COUNT;
//			stateLikelihood.calculateStatesLogLikelihood(spectra, kOffset, allStateLogLikelihood);
//			stateLikelihood.calculateStoredStatesLogLikelihood(spectra, kOffset, allStoredStateLogLikelihood);
//		}
//		for (int i : mapToSrp) {
//			int state = allSrpState2D[i][k];
//			for (int j : allSpectrumIndexs) {
////				currentLogLikelihood = updateLikelihoodAtIJK(i, j, state, allStateLogLikelihood2D[j],
////						allStoredStateLogLikelihood2D[j], currentLogLikelihood);
//				if (state < STATE_COUNT) {
//					currentLogLikelihood = updateLikelihoodAtIJK(i, j, j*STATE_COUNT+state,
//							currentLogLikelihood);
//				}
//			}
//		}
		char[] oldChars = new char[allSpectrumIndexs.length]; 
		char[] newChars = new char[allSpectrumIndexs.length];
		int[] shortList = new int[allSpectrumIndexs.length];
		int count = 0;
		for (int j : allSpectrumIndexs) {
			Haplotype haplotype = alignmentModel.getHaplotype(j);
			oldChars[j] = haplotype.getStoredChar(k);
			newChars[j] = haplotype.getChar(k);
			if(oldChars[j] != newChars[j]){
				shortList[count++] = j;
			}
		}
//		shortList.toArray(a)

		for (int i : mapToSrpArray[k]){
			char srpChar = allSrpChar2D[i][k];
			for (int s = 0; s < count; s++) {
				int j = shortList[s];
				int deltaDist = calculateDeltaDist(srpChar, newChars[j], oldChars[j]);
				allDists[i][j] += deltaDist;
			}
//		}
//		for (int i : mapToSrpArray[k]){
			currentLogLikelihood = updateEachSrpAtI(i, currentLogLikelihood);
		}
		return currentLogLikelihood;
	}
	
	protected double calculateSrpLikelihoodSubColumn() {

//		int k = operationRecord.getSingleIndex();
//		int[] allSpectrumIndexs = operationRecord.getAllSpectrumIndexs();
		
		double currentLogLikelihood = getStoredLogLikelihood();

		return currentLogLikelihood;
	}

	

	protected double calculateSrpLikelihoodRecombination() {
//		System.out.println("Cal srp recombination");

		int[] twoSpectrums = operationRecord.getRecombinationSpectrumIndex();
		int[] twoPositions = operationRecord.getRecombinationPositionIndex();

//		Spectrum[] spectrums = new Spectrum[] {
//				spectrumModel.getHaplotype(twoSpectrums[0]),
//				spectrumModel.getHaplotype(twoSpectrums[1]) };
	
		int j0 = twoSpectrums[0];
		int j1 = twoSpectrums[1];
		int length = twoPositions[1] - twoPositions[0];
		
		Haplotype haplotype0 = alignmentModel.getHaplotype(j0); 
		Haplotype haplotype1 = alignmentModel.getHaplotype(j1);
		
		int[] siteIndexs = new int[length];
		for (int k = twoPositions[0], s=0; k < twoPositions[1]; k++, s++) {
			siteIndexs[s] = k;
//			ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);
//			for (int i : mapToSrp) {
//				srpSwitch[i] = true;
//			}
		}

//		for (int k : siteIndexs) {
//			
//			SpectraParameter spectra0 = spectrum0.getSpectra(k);
//			SpectraParameter spectra1 = spectrum1.getSpectra(k);
//			int kOffset = k*STATE_COUNT;
//			stateLikelihood.calculateStoredStatesLogLikelihood(spectra0, kOffset, allStateLogLikelihood);
//			stateLikelihood.calculateStoredStatesLogLikelihood(spectra1, kOffset, allStoredStateLogLikelihood);
//			
////			stateLikelihood.calculateStoredStatesLogLikelihood2D(spectra0, allStateLogLikelihood2D[s]);
////			stateLikelihood.calculateStoredStatesLogLikelihood2D(spectra1, allStoredStateLogLikelihood2D[s]);
//			spectra0.setStateLikelihood(allStoredStateLogLikelihood, kOffset);
//			spectra1.setStateLikelihood(allStateLogLikelihood, kOffset);
//		
//		}

		int multihere;
		double currentLogLikelihood = getStoredLogLikelihood();
		if(multiType == MultiType.BitSet){
			recalculateBitSet(siteIndexs);
			
			int count = 0;
			for (int i = bitSet.nextSetBit(0); i >= 0; i = bitSet.nextSetBit(i+1)) {
				srpIndex[count++] = i;
//				count++;
				currentLogLikelihood = updateLikelihoodRecomb(i, j0, j1, siteIndexs, haplotype0, haplotype1, currentLogLikelihood);
//				currentLogLikelihood = updateLikelihoodAtIJ(i, j0,
//						siteIndexs, haplotype0, currentLogLikelihood);
//				currentLogLikelihood = updateLikelihoodAtIJ(i, j1,
//						siteIndexs, haplotype1, currentLogLikelihood);
			}
			srpIndexCount = count;

			
		}
		else if(multiType==MultiType.Array){

			for (int s : siteIndexs) {
				for (int i : mapToSrpArray[s]){
					srpSwitch[i] = true;
				}
			}

			for (int i = 0; i < srpSwitch.length; i++) {
				if(srpSwitch[i]){

					
				}
			}
	
		}
		

		
	
//		totalLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return currentLogLikelihood;

	}






	private double updateLikelihoodAtIJ(int i, int j, int[] siteIndexs, 
			Haplotype haplotype, double currentLogLikelihood) {

		updateAllDists(i, j, siteIndexs, haplotype);
		currentLogLikelihood = updateEachSrpAtI(i, currentLogLikelihood);

		return currentLogLikelihood;
	}
	
	int deltaCount;
	private double updateLikelihoodRecomb(int i, int j0, int j1, int[] siteIndexs, 
			Haplotype haplotype0, Haplotype haplotype1, double currentLogLikelihood) {
		deltaCount=0;
		updateAllDists(i, j0, siteIndexs, haplotype0);
//		currentLogLikelihood = updateEachSrpAtI(i, currentLogLikelihood);
		updateAllDists(i, j1, siteIndexs, haplotype1);
		currentLogLikelihood = updateEachSrpAtI(i, currentLogLikelihood);
		
//		updateLikelihoodAtIJ(i, j0, siteIndexs, haplotype0, currentLogLikelihood)
//		System.out.println(deltaCount +"\t"+ i);
		return currentLogLikelihood;
		
	}



	private double updateEachSrpAtI(int i, double currentLogLikelihood, int... dists) {
		
		
		int srpLength = srpMap.getSrpLength(i);
		double[] logPD = scaledLogBinomialDesnity.get(srpLength);
		
		double srpSum = 0;//~60 vs 90 for multi
//		liS.reset();//LiS
		try{
		for (int hj = 0; hj < sequenceCount; hj++) {
			srpSum += logPD[allDists[i][hj]];
//			liS.add(logPD[allDists[i][hj]]);//LiS
//			srpSum += logPD[hj];
		}
		}
		catch (ArrayIndexOutOfBoundsException e){
			System.out.println(i);
			System.out.println(Arrays.toString(allDists[i]));
			for (int hj = 0; hj < sequenceCount; hj++) {
				System.out.println(hj);
				System.out.println(logPD[allDists[i][hj]]);
				srpSum += logPD[allDists[i][hj]];
			}
		}
//		BigDecimal bd = new BigDecimal(srpSum);
//		bd.
//		int oldDist = dists[0];
//		int newDist = dists[1];
//		sumScaledSrpLogLikelihood[i] = storedSumScaledSrpLogLikelihood[i] - logPD[oldDist]+ logPD[newDist];
//		srpSum = sumScaledSrpLogLikelihood[i];
		
		currentLogLikelihood -= eachSrpLogLikelihood[i];
		eachSrpLogLikelihood[i] = LikelihoodScaler.getLogLikelihood(srpSum, LOG_C);
//		eachSrpLogLikelihood[i] = liS.getLogLikelihood();//LiS
		currentLogLikelihood += eachSrpLogLikelihood[i];
		return currentLogLikelihood;
		
	}



	private void updateAllDists(int i, int j, int[] siteIndexs, Haplotype haplotype) {
		
		char[] srpChars = allSrpChar2D[i];
		int deltaDist = 0;
		for (int k : siteIndexs) {
			char oldChar = haplotype.getStoredChar(k);
			char newChar = haplotype.getChar(k);
			if(newChar!= oldChar){ 
				deltaDist += calculateDeltaDist(srpChars[k], newChar, oldChar);
				deltaCount++;
			}
		}
//		int oldDist = storedAllDists[i][j];
//		int newDist = storedAllDists[i][j] + deltaDist;
		allDists[i][j] += deltaDist;
		if(allDists[i][j]<0){
			System.out.println(allDists[i][j] -= deltaDist);
			System.out.println(allDists[i][j] += deltaDist);
			System.out.println(Arrays.toString(siteIndexs));
			for (int k : siteIndexs) {
				char oldChar = haplotype.getStoredChar(k);
				char newChar = haplotype.getChar(k);
				System.out.println("site:"+k +"\t"+ i +"\t"+ j);
				System.out.println(srpChars[k] +"\t"+ newChar +"\t"+ oldChar);
				System.out.println(calculateDeltaDist(srpChars[k], newChar, oldChar) +"\t"+ allDists[i][j]);
			}
//			for (int j = 0; j < sequenceCount; j++) {
				
				String srp = srpMap.getSrpFragment(i);//srpArray[i]
				int start = srpMap.getSrpStart(i);
				int end = srpMap.getSrpEnd(i);
				
				int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(j).getSequenceString());
				System.out.println(start +"\t"+ end +"\t"+ srp);
				System.out.println(alignmentModel.getHaplotype(j).getSequenceString());
				System.out.println("trueDist:\t"+dist);
//				liS.scaleLogProb(logPD[dist]);
//			}	

		}
		
	}



	@Override
	protected void storeState() {

//		System.arraycopy(eachSrpLogLikelihood, 0, storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		storedLogLikelihood = logLikelihood;
//		storedBigDecLikelihood = bigDecLikelihood;
//		storeEverything();
		
		OperationType operation = operationRecord.getOperation();
//		int spectrumIndex;
//		int siteIndex; = spectrumOperationRecord.getAllSiteIndexs()[0];

		int j;
		int k;
		if(DEBUG){
			System.out.println("StoreState in ShortReadsSpectrumLikelihood:\t"+operation);
		}
		switch (operation) {
		case NONE:
			break;
		case FULL:
			storeEverything();
			break;

		case COLUMN:
		
		case SWAP_SUBCOLUMN:
			k= operationRecord.getSingleIndex();
			
			for (int i : mapToSrpArray[k]){
				storeI(i);
				for (j = 0; j < sequenceCount; j++) {
					storeIJ(i, j);
				}
			}

			break;
		case SINGLE:
		
			j = operationRecord.getSpectrumIndex();
			k = operationRecord.getSingleIndex();//AllSiteIndexs()[0];
			
			for (int i : mapToSrpArray[k]){
				storeI(i);
				storeIJ(i, j);
				
//				storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
//				storedScaledSpectrumLogLikelihood[i][j] = scaledSpectrumLogLikelihood[i][j];
//				storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
//				storedSumSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i];

			}
			break;
			
		case MULTI:
		
			j = operationRecord.getSpectrumIndex();
			int[] siteIndexs = operationRecord.getAllSiteIndexs();
			int storeMulti;
			if(multiType==MultiType.BitSet){
				for (int s = 0; s < srpIndexCount; s++) {
					int i = srpIndex[s];
					storeI(i);
					storeIJ(i, j);
				}
			}
			
			else if(multiType==MultiType.Array){
				for (int i = 0; i < srpSwitch.length; i++) {
					if (srpSwitch[i]) {
						storeI(i);
						storeIJ(i, j);
						srpSwitch[i] = false;
					}
				}
			}
			else if(multiType==MultiType.Hash){
				for (int i : allSrpPos) {
					storeI(i);
					storeIJ(i, j);
				}
			}
			else if(multiType==MultiType.All){
				for (int kk : siteIndexs) {
					for (int i : mapToSrpArray[kk]){
						storeI(i);
						storeIJ(i, j);
					}
				}
			}
						
			break;
		case RECOMBINATION:
			int[] twoSpectrums = operationRecord.getRecombinationSpectrumIndex();
//			int[] twoPositions = spectrumOperationRecord.getRecombinationPositionIndex();

//			for (int i = 0; i < srpSwitch.length; i++) {
//				if (srpSwitch[i]) {
//					storeI(i);
//					storeIJ(i, twoSpectrums[0]);
//					storeIJ(i, twoSpectrums[1]);
//					srpSwitch[i] = false;
//				}
//			}
			
			if(multiType==MultiType.BitSet){
				for (int s = 0; s < srpIndexCount; s++) {
					int i = srpIndex[s];
					storeI(i);
					storeIJ(i, twoSpectrums[0]);
					storeIJ(i, twoSpectrums[1]);

				}
			}
			
			else if(multiType==MultiType.Array){
				for (int i = 0; i < srpSwitch.length; i++) {
					if (srpSwitch[i]) {
						storeI(i);
						storeIJ(i, twoSpectrums[0]);
						storeIJ(i, twoSpectrums[1]);
						srpSwitch[i] = false;
					}
				}
			}
			
			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: "+operation +"\tin"+ShortReadsHaplotypeLikelihood.class.getSimpleName() );
//			storeEverything();
//			break;
			
		}


//		long time2 = System.currentTimeMillis();
//		time += (time2-time1);
		
		
//		System.out.println(spectrumIndex +"\t"+ siteIndex);
//		ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(siteIndex);
//		for (int i : mapToSrp) {
//			storedAllLogLikelihood[i][spectrumIndex][siteIndex] = allLogLikelihood[i][spectrumIndex][siteIndex];
//			storedSpectrumLogLikelihood[i][spectrumIndex] = spectrumLogLikelihood[i][spectrumIndex];
//		}
		
//		for (int i = 0; i < srpCount; i++) {
//			System.err.println(i +"\t"+ allLogLikelihood[i][spectrumIndex][siteIndex] +"\t"+ storedAllLogLikelihood[i][spectrumIndex][siteIndex]);
//			storedAllLogLikelihood[i][spectrumIndex][siteIndex] = allLogLikelihood[i][spectrumIndex][siteIndex];
//			for (int j = 0; j < spectrumCount; j++) {
//				for (int k = 0; k < spectrumLength; k++) {
//					if(allLogLikelihood[i][j][k] != storedAllLogLikelihood[i][j][k]){
//						System.out.println("DIFFLI:"+i +" "+j+" "+" "+k+
//								" "+ allLogLikelihood[i][j][k] +" "+
//								storedAllLogLikelihood[i][j][k]);
//					}
//				}
//				System.arraycopy(allLogLikelihood[i][j], 0, storedAllLogLikelihood[i][j], 0, spectrumLength);
//			}
//			storedAllLogLikelihood[i][spectrumIndex][siteIndex] = allLogLikelihood[i][spectrumIndex][siteIndex];
//			
//		}
		
		
//		System.err.println("SR likelihood store: " + logLikelihood +"\t"+ storedLogLikelihood);
//		System.err.println(allLogLikelihood);
	}
//	private void storeMultiDelta(SpectrumOperationRecord spectrumOperationRecord) {
//		
//		
//	}


	private void storeI(int i) {
		storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
		storedSumScaledSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i];
//		storedSumBigDecSrp[i] = sumBigDecSrp[i];
//		storedBigDecEachSrpLogLikelihood[i] = bigDecEachSrpLogLikelihood[i];
	}


	private void storeIJ(int i, int j) {
		storedAllDists[i][j] = allDists[i][j];
		
	}


	//	private void storeMultiDelta(SpectrumOperationRecord spectrumOperationRecord) {
	//		
	//		
	//	}
	
	
	private void storeEverything() {

		System.arraycopy(eachSrpLogLikelihood, 0, storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		System.arraycopy(sumScaledSrpLogLikelihood, 0, storedSumScaledSrpLogLikelihood, 0, sumScaledSrpLogLikelihood.length);
		
//		System.arraycopy(sumBigDecSrp, 0, storedSumBigDecSrp, 0, sumBigDecSrp.length);
//		System.arraycopy(bigDecEachSrpLogLikelihood, 0, storedBigDecEachSrpLogLikelihood, 0, bigDecEachSrpLogLikelihood.length);
		
		for (int i = 0; i < allDists.length; i++) {
			System.arraycopy(allDists[i], 0, storedAllDists[i], 0, sequenceCount);
		}
		
	}


	@Override
	protected void restoreState() {
//		long time1 = System.currentTimeMillis();
		
//		System.err.println("SR likelihood restore: " + logLikelihood +"\t"+ storedLogLikelihood);
		logLikelihood = storedLogLikelihood;
		
//		BigDecimal temp2 = storedBigDecLikelihood;
//		storedBigDecLikelihood = bigDecLikelihood;
//		bigDecLikelihood = temp2;
		
//		restoreEverything();
//		System.arraycopy(storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		
//		OperationRecord spectrumOperationRecord = alignmentModel.getOperationRecord();
		OperationType operation = operationRecord.getOperation();
//		int spectrumIndex;
//		int siteIndex = spectrumOperationRecord.getAllSiteIndexs()[0];
//		ArrayList<Integer> mapToSrp;
//		int[] siteIndexs;
		int j;
		int k;
		if(DEBUG){
			System.out.println("RestoreState in ShortReadsSpectrumLikelihood:\t"+operation);
		}
		switch (operation) {
		case NONE:
			
			break;
		case FULL:
//			System.out.println("Ever call this? restoreEverything??");
			restoreEverything();
			break;
		case COLUMN:
		
			k = operationRecord.getSingleIndex();
			for (int i : mapToSrpArray[k]){
				restoreI(i);
				for (j = 0; j < sequenceCount; j++) {
					restoreIJ(i, j);
				}
			}

			break;
		case SINGLE:

			j = operationRecord.getSpectrumIndex();
			k = operationRecord.getSingleIndex();//AllSiteIndexs()[0];
//			mapToSrp = srpMap.getMapToSrp(k);
//			spectrumModel.getSpectrum(j).getSpectra(k).restoreState();
//			for (int i : mapToSrp) {
			for (int i : mapToSrpArray[k]){
				restoreI(i);
				restoreIJ(i, j);
			}
			
			
		case MULTI:
		
			j = operationRecord.getSpectrumIndex();
			int[] siteIndexs = operationRecord.getAllSiteIndexs();
			int restoreMulti;
			if(multiType==MultiType.BitSet){

				for (int s = 0; s < srpIndexCount; s++) {
					int i = srpIndex[s];
					restoreI(i);
					restoreIJ(i, j);
				}
			}

			else if(multiType==MultiType.Array){

				for (int i = 0; i < srpSwitch.length; i++) {
					if (srpSwitch[i]) {
						restoreI(i);
						restoreIJ(i, j);
					}
				}
			}

			else if(multiType==MultiType.Hash){
				for (int i : allSrpPos) {
					restoreI(i);
					restoreIJ(i, j);
				}
				
			}
			else if(multiType==MultiType.All){
				for (int kk : siteIndexs) {
					for (int i : mapToSrpArray[kk]){
						restoreI(i);
						restoreIJ(i, j);
					}
				}

			}
			
			break;
		case RECOMBINATION:
			int[] twoSpectrumsIndex = operationRecord.getRecombinationSpectrumIndex();

			if(multiType==MultiType.BitSet){
				for (int s = 0; s < srpIndexCount; s++) {
					int i = srpIndex[s];
					restoreI(i);
					restoreIJ(i, twoSpectrumsIndex[0]);
					restoreIJ(i, twoSpectrumsIndex[1]);
				}
			}

			else if(multiType==MultiType.Array){

				for (int i = 0; i < srpSwitch.length; i++) {
					if (srpSwitch[i]) {
						restoreI(i);
						restoreIJ(i, twoSpectrumsIndex[0]);
						restoreIJ(i, twoSpectrumsIndex[1]);
					}
				}
			}

			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: "+operation +"\tin"+ShortReadsHaplotypeLikelihood.class.getSimpleName() );

		}
//		long time2 = System.currentTimeMillis();
//		time += (time2-time1);
//
//		
//		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
//		int spectrumIndex = record.getSpectrumIndex();
//		int siteIndex = record.getSiteIndex();
//
//		ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(siteIndex);
//		for (int i : mapToSrp) {
//			allLogLikelihood[i][spectrumIndex][siteIndex] = storedAllLogLikelihood[i][spectrumIndex][siteIndex];
//			spectrumLogLikelihood[i][spectrumIndex] = storedSpectrumLogLikelihood[i][spectrumIndex];
//		}
//		
		
//		for (int i = 0; i < srpCount; i++) {
//			System.err.println(i +"\t"+ allLogLikelihood[i][spectrumIndex][siteIndex] +"\t"+ storedAllLogLikelihood[i][spectrumIndex][siteIndex]);
//			allLogLikelihood[i][spectrumIndex][siteIndex] = storedAllLogLikelihood[i][spectrumIndex][siteIndex];
//			for (int j = 0; j < spectrumCount; j++) {
//				
//				System.arraycopy(storedAllLogLikelihood[i][j], 0, allLogLikelihood[i][j], 0, spectrumLength);
//			}
//			
//		}
//		System.err.println(
//				spectrumIndex +"\t"+ siteIndex +"\t"
//						+Arrays.toString(spectrumModel.getSpectrum(spectrumIndex)
//						.getFrequencies(siteIndex)));
//
//		System.err.println("SR likelihood restore: " + logLikelihood +"\t"+ storedLogLikelihood);
		
	}


	private void restoreI(int i) {
		eachSrpLogLikelihood[i] = storedEachSrpLogLikelihood[i];
		sumScaledSrpLogLikelihood[i] = storedSumScaledSrpLogLikelihood[i]; 
		
//		temp = storedstateLikelihood;
//    	storedstateLikelihood = stateLikelihood;
//    	stateLikelihood = temp;
//		BigDecimal temp = storedSumBigDecSrp[i];
//		storedSumBigDecSrp[i] = sumBigDecSrp[i]; 
//		sumBigDecSrp[i] = temp;
//		
//		BigDecimal temp2 = storedBigDecEachSrpLogLikelihood[i];
//		storedBigDecEachSrpLogLikelihood[i] = bigDecEachSrpLogLikelihood[i];
//		bigDecEachSrpLogLikelihood[i] = temp2;
	}


	private void restoreIJ(int i, int j) {
		allDists[i][j] = storedAllDists[i][j];
		
	}


	private void restoreEverything() {
		
		System.arraycopy(storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		System.arraycopy(storedSumScaledSrpLogLikelihood, 0, sumScaledSrpLogLikelihood, 0, sumScaledSrpLogLikelihood.length);
//		System.arraycopy(storedSumBigDecSrp, 0, sumBigDecSrp, 0, sumBigDecSrp.length);
//		System.arraycopy(storedBigDecEachSrpLogLikelihood, 0, bigDecEachSrpLogLikelihood, 0, bigDecEachSrpLogLikelihood.length);
		for (int i = 0; i < allDists.length; i++) {
			System.arraycopy(storedAllDists[i], 0, allDists[i], 0, sequenceCount);
		}
		
	}


	public double[] unittestMethodGetEachLikelihood() {
		double[] copyOfValues = new double[eachSrpLogLikelihood.length];
		System.arraycopy(eachSrpLogLikelihood, 0, copyOfValues, 0,
				copyOfValues.length);
		return copyOfValues;
	}
	
	public double calculateSrpLikelihoodFullMaster() {


//		System.out.println("calculateSrpLikelihoodMaster");
		double logLikelihood = 0;
		for (int i = 0; i < srpCount; i++) {
			
			String srp = srpMap.getSrpFragment(i);//srpArray[i]
			int start = srpMap.getSrpStart(i);
			int end = srpMap.getSrpEnd(i);
			
			double[] logPD = scaledLogBinomialDesnity.get(srpMap.getSrpLength(i));

			liS.reset();
			for (int j = 0; j < sequenceCount; j++) {

				int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(j).getSequenceString());
				liS.add(logPD[dist]);
//				liS.scaleLogProb(logPD[dist]);
			}	
			
			double eachSrp = liS.getLogLikelihood();
			if(DEBUG){
				if(eachSrp != eachSrpLogLikelihood[i]){
					System.out.println(i +"\t"+ eachSrp +"\t"+ eachSrpLogLikelihood[i]);
				}
			}
			logLikelihood += eachSrp;

		}
//		 = StatUtils.sum(eachSrpLogLikelihood);
		
//for (int i = 0; i < srpCount; i++) {
//	System.out.println(Arrays.toString(allDists[i]));
//}
//System.out.println("==");

//		double logLikelihood2 = calculateSrpLikelihoodFull();
//		if(logLikelihood != logLikelihood2){
//			System.out.println(logLikelihood +"\t"+ logLikelihood2 +"\t"+ (logLikelihood-logLikelihood2));
//		}
//		
		
		return logLikelihood;
	}

	@Override
	public void makeDirty() {
		alignmentModel.resetOperation();
		likelihoodKnown = false;
	
	}


//
//	protected double calculateSrpLikelihoodFullBig() {
//	
//	//		System.out.println("calculateSrpLikelihood_Full");
//			double logLikelihood = 0;
//			bigDecLikelihood = new BigDecimal("0");
//			for (int i = 0; i < srpCount; i++) {
//				
//				String srp = srpMap.getSrpFragment(i);//srpArray[i]
//				int start = srpMap.getSrpStart(i);
//				int end = srpMap.getSrpEnd(i);
//				
//				double[] logPD = scaledLogBinomialDesnity.get(srpMap.getSrpLength(i));
//				BigDecimal[] bigDec = bigDecBinomialDesnity.get(srpMap.getSrpLength(i));
//				liS.reset();
//				sumScaledSrpLogLikelihood[i] = 0;
//				sumBigDecSrp[i] = new BigDecimal("0");
//				for (int j = 0; j < sequenceCount; j++) {
//	
//					int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(j));
//					allDists[i][j]=dist;
//					liS.add(logPD[dist]);
////					sumBigDecSrp[i] = sumBigDecSrp[i].add(bigDec[dist]);
//	//				System.out.println(sumBigDecSrp[i]);
//	//				sumScaledSrpLogLikelihood[i] += logPD[allDists[i][j]];
//	//				storedAllLog2D[i][j] = logPD[dist];
//	//				liS.scaleLogProb(logPD[dist]);
//				}	
//				sumScaledSrpLogLikelihood[i] = liS.getSumScaledLikelihood();
//				eachSrpLogLikelihood[i] = liS.getLogLikelihood() ;
//				
//				
//				bigDecEachSrpLogLikelihood[i] = BigFunctions.ln(sumBigDecSrp[i], 300);
//				eachSrpLogLikelihood[i] = bigDecEachSrpLogLikelihood[i].doubleValue();
//	//			System.out.println(eachSrpLogLikelihood[i] +"\t"+ sumBigDecSrp[i].toPlainString());
//	//			System.out.println("inFull: "+sumScaledSrpLogLikelihood[i] +"\t"+ eachSrpLogLikelihood[i]);
//				bigDecLikelihood = bigDecLikelihood.add(bigDecEachSrpLogLikelihood[i]);
//				logLikelihood += eachSrpLogLikelihood[i];
//			}
//			System.out.println(logLikelihood +"\t"+ bigDecLikelihood.toPlainString());
//			return logLikelihood;
//		}


//
//	protected double calculateSrpLikelihoodSingleBig() {
//	
//			int j = operationRecord.getSpectrumIndex(); 
//			int k = operationRecord.getSingleIndex();//AllSiteIndexs()[0];
//			
//			double currentLogLikelihood = getStoredLogLikelihood();
//	
//			Haplotype haplotype = alignmentModel.getHaplotype(j);
//			char oldChar = haplotype.getStoredChar(k);
//			char newChar = haplotype.getChar(k);
//			
//	//		System.out.print("\t"+ currentLogLikelihood +"\t" );
//			if(newChar!= oldChar){ 
//				for (int i : mapToSrpArray[k]){
//					char srpChar = allSrpChar2D[i][k];
//					int deltaDist = calculateDeltaDist(srpChar, newChar, oldChar);
//					
//					if (deltaDist != 0) {
//						int srpLength = srpMap.getSrpLength(i);
//						BigDecimal[] tempBD = bigDecBinomialDesnity.get(srpLength);
//						double[] logPD = scaledLogBinomialDesnity.get(srpLength);
//						
//						double sum = 0;
//						
//	//					double sum = 0;
//	//					BigDecimal bd = new BigDecimal("0");
//	//					bd.setScale(100);
//						for (int l = 0; l < sequenceCount; l++) {
//							sum += logPD[allDists[i][l]];
//	//						bd = bd.add(BigDecimal.valueOf( logPD[allDists[i][l]] ));
//						}
//						if(deltaDist == -1){
//							System.out.println(deltaDist +"\t"+  srpChar +"\t"+ newChar +"\t"+ oldChar +"\t"+ cheat[j].charAt(k));
//							System.out.println("Old\t"+currentLogLikelihood +"\t"+ eachSrpLogLikelihood[i] +"\t"+ sumBigDecSrp[i].doubleValue() +"\t"+  sumBigDecSrp[i].toPlainString());
//						}
//						sumBigDecSrp[i] = sumBigDecSrp[i].subtract(tempBD[allDists[i][j]]);
//						
//						allDists[i][j] += deltaDist;
//						
//						sumBigDecSrp[i] = sumBigDecSrp[i].add(tempBD[allDists[i][j]]);
//						sum = 0;
//						for (int l = 0; l < sequenceCount; l++) {
//							sum += logPD[allDists[i][l]];
//	//						bd = bd.add(BigDecimal.valueOf( logPD[allDists[i][l]] ));
//						}
//						bigDecLikelihood = bigDecLikelihood.subtract(bigDecEachSrpLogLikelihood[i]);
//						currentLogLikelihood -= eachSrpLogLikelihood[i];
//	//					BigDecimal bd = BigFunctions.ln(sumBigDecSrp[i], sumBigDecSrp[i].scale());
//						
//						
//						bigDecEachSrpLogLikelihood[i] = BigFunctions.ln(sumBigDecSrp[i], 300);
//						eachSrpLogLikelihood[i] = bigDecEachSrpLogLikelihood[i].doubleValue();
//						bigDecLikelihood = bigDecLikelihood.add(bigDecEachSrpLogLikelihood[i]);
//						
//						
//	//					eachSrpLogLikelihood[i] = bd.doubleValue();
//	//					eachSrpLogLikelihood[i] = liS.getLogLikelihood();//LiS
//						currentLogLikelihood += eachSrpLogLikelihood[i];
//						
//	//					currentLogLikelihood = updateEachSrpAtI(i, currentLogLikelihood, (allDists[i][j]-deltaDist), (allDists[i][j]));
//						if(deltaDist == -1){
//	//						double newSum = 0;
//	//						double allLogSum = 0;
//	//						BigDecimal newBd = new BigDecimal("0");
//	//						newBd.setScale(100);
//	////						int srpLength = srpMap.getSrpLength(i);
//	////						double[] logPD = scaledLogBinomialDesnity.get(srpLength);
//	//						for (int l = 0; l < sequenceCount; l++) {
//	////							newSum += allDists[i][l];
//	//							newSum += logPD[allDists[i][l]];
//	//							BigDecimal tempB = BigDecimal.valueOf( logPD[allDists[i][l]] );
//	//							
//	//							newBd = newBd.add(tempB);
//	//							
//	//						}
//							System.out.println("new\t"+currentLogLikelihood +"\t"+ eachSrpLogLikelihood[i]  +"\t"+ sumBigDecSrp[i].doubleValue() +"\t"+  sumBigDecSrp[i].toPlainString());
//							System.out.println(LikelihoodScaler.getLogLikelihood(sum, LOG_C) +"\t"+ bigDecEachSrpLogLikelihood[i].doubleValue() +"\t"+ bigDecEachSrpLogLikelihood[i].toPlainString());
//							System.out.println(currentLogLikelihood +"\t"+ bigDecLikelihood.doubleValue() +"\t"+ bigDecLikelihood.toPlainString());
//	//						System.out.println("===\t"+ newBd.doubleValue());
//	//						double log2 = newSum*LOG_ERROR_RATE+(srpLength*sequenceCount-newSum)*LOG_ONE_MINUS_ERROR_RATE;
//	//						System.out.println(sum +"\t"+ bd.toString() +"\n"+ newSum +"\t"+ newBd.toString() );  
//	//						System.out.println(liS.getLogLikelihood(bd.doubleValue(), LOG_C) +"\t"+ liS.getLogLikelihood(newBd.doubleValue(), LOG_C) );
//	//						System.out.println(liS.getLogLikelihood(Double.valueOf(bd.toPlainString()), LOG_C) +"\t"+ liS.getLogLikelihood(Double.valueOf(newBd.toPlainString()), LOG_C) );
//	//						System.out.println(BigFunctions.ln(bd, 100).toPlainString() +"\n"+ BigFunctions.ln(newBd, 100).toPlainString());
//	//						System.out.println((BigFunctions.ln(bd, 100).doubleValue())-LOG_C +"\n"+ (BigFunctions.ln(newBd, 100).doubleValue()-LOG_C));
//	//						System.out.println(allDists[i][j] +"\t"+  logPD[allDists[i][j]] +"\t"+ logPD[ (allDists[i][j]-deltaDist  )]);
//							System.out.println(sumScaledSrpLogLikelihood[i] +"\t"+  storedSumScaledSrpLogLikelihood[i] +"\t"+ storedBigDecLikelihood.toPlainString());
//							
//							System.out.println();
//	//						System.out.println(LOG_C);
//						}
//					}
//	
//				}
//			}
//			currentLogLikelihood = bigDecLikelihood.doubleValue();
//	//		System.out.println(currentLogLikelihood);
//			return currentLogLikelihood;
//		}

//
//		double logLikelihood = 0;
//		double spectrumLogLikelihood = 0;
//		double stateLogLikelihood = 0;
//		
//		for (int i = 0; i < srpCount; i++) {
//
//			String fullSrp = srpMap.getSrpFull(i);
//			int start = srpMap.getSrpStart(i);
//			int end = srpMap.getSrpEnd(i);
//			
//			liS.reset();
//			for (int j = 0; j < spectrumCount; j++) {
//
//				Haplotype haplotype = haplotypeModel.getHaplotype(j);
//				spectrumLogLikelihood = 0;
//				for (int k = start; k < end; k++) {
//					double[] frequencies = spectrum.getFrequenciesAt(k);
//					int state = getStateAtK(fullSrp, k);
//					if(state<STATE_COUNT){
////						stateLogLikelihood = caluclateStateLogLikelihood(frequencies[state]);
//						stateLogLikelihood = stateLikelihood.caluclateStateLogLikelihood(frequencies[state]);
////						double likelihood = frequencies[state] * NOT_ERROR_RATE
////								+ (1 - frequencies[state]) * ERROR_RATE;
////						stateLogLikelihood = Math.log(likelihood);
//					}
//					else{
//						stateLogLikelihood = MIN_LOG_LIKELIHOOD;
//					}
//					
//					spectrumLogLikelihood += stateLogLikelihood;
//					
//				}
//				liS.addLogProb(spectrumLogLikelihood);
//				
//
//			}	
//			logLikelihood += liS.getLogLikelihood();
//		}
//
////		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
//		if(DEBUG){
//			if(logLikelihood != this.logLikelihood){
//				System.out.println(logLikelihood +"\t"+ this.logLikelihood +"\t"+ getStoredLogLikelihood());
//			
//			
//			
//				OperationRecord record = haplotypeModel.getOperationRecord();
////				int[] siteIndexs = record.getAllSiteIndexs();
//				int j= record.getSpectrumIndex(); 
//				Haplotype haplotype = haplotypeModel.getHaplotype(j);
//				
////				stateLikelihood.calculateStatesLogLikelihood(spectra, allStateLogLikelihood[k]);
////				stateLikelihood.calculateStoredStatesLogLikelihood(spectra, allStoredStateLogLikelihood[k]);
//				
//				for (int s = 0; s < spectrumLength; s++) {
//					int k = s;
//					SpectraParameter spectra = spectrum.getSpectra(k);
////					int kOffset = k*STATE_COUNT;
//					double[] stateLn = spectra.getStateLikelihood();
////					stateLikelihood.calculateStatesLogLikelihood(spectra, 0, stateLn);
////					stateLikelihood.calculateStoredStatesLogLikelihood(spectra, kOffset, allStoredStateLogLikelihood);
////					System.out.println(Arrays.toString());
//					
//					double[] frequencies = spectrum.getFrequenciesAt(k);
//					double[] stateLn2 = new double[4];
//					for (int state = 0; state < 4; state++) {
//						stateLn2[state] = stateLikelihood.caluclateStateLogLikelihood(frequencies[state]);
//						if(stateLn2[state] != stateLn[state]){
//							System.out.println(s);
//							System.out.println(Arrays.toString(stateLn));
//							System.out.println(Arrays.toString(stateLn2));
//						}
//					}
////					System.out.println(Arrays.toString(stateLn2));	
//				}
////				System.exit(-1);
//			
//			
//			
//			
//			
//			
//			
//			}
//		}
//		return logLikelihood;
//	}
//

}
