package srp.spectrum.likelihood;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.ArithmeticUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.Operation;
import srp.haplotypes.ShortRead;
import srp.haplotypes.SwapInfo;
import srp.likelihood.LikelihoodScaler;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import srp.spectrum.SpectrumOperationRecord;
import dr.evolution.datatype.DataType;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

public class ShortReadsSpectrumLikelihood  extends AbstractModelLikelihood {

	
    public static final String SHORT_READ_LIKELIHOOD = "ShortReadSpectrumLikelihood";
	public static final String NAME = SHORT_READ_LIKELIHOOD;
	public static final double ERROR_RATE = 0.0107;
	public static final double NOT_ERROR_RATE = 1-ERROR_RATE;
	public static final double LOG_ERROR_RATE = Math.log(ERROR_RATE);
//	public static final double LOG_ONE_MINUS_ERROR_RATE = Math.log(1-ERROR_RATE);
	public static final double C = 1e-200;
	public static final double LOG_C = Math.log(C);
	
	private static final double EVALUATION_TEST_THRESHOLD = 1e-8;
	private static final int TWO = 2;
	
	//	public static final int[] NULL_SWAPINFO = HaplotypeModel.NULL_SWAPINFO;
	protected boolean likelihoodKnown;
	
	private int spectrumLength;
	private int spectrumCount;
	private int srpCount;

	private double logLikelihood;
	private double storedLogLikelihood;

	private double[] eachSrpLikelihood;
	private double[] storedEachSrpLikelihood;

	private double[][][] allLogLikelihood;
	private double[][][] storedAllLogLikelihood;

//	private int[] allDelta;


//	private HashMap<Integer, double[]> logBinomialDesnity;
//	private HashMap<Integer, double[]> scaledLogBinomialDesnity;

	private AlignmentMapping aMap;
	
//	private HaplotypeModel haplotypeModel;
	private LikelihoodScaler liS;
//	private SwapInfo swapInfo;
	

	private SpectrumOperation operation;
	private SpectrumAlignmentModel spectrumModel;
	private DataType dataType;
	private int stateCount;
	private double[][] spectrumLogLikelihood;
	private double[][] storedSpectrumLogLikelihood;
	
	
	@Override
	public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

    public ShortReadsSpectrumLikelihood(String name) {

        super(SHORT_READ_LIKELIHOOD);
        likelihoodKnown = false;
        setId(SHORT_READ_LIKELIHOOD);
        

    }


	public ShortReadsSpectrumLikelihood(SpectrumAlignmentModel spectrumModel){
		this(SHORT_READ_LIKELIHOOD);
		this.spectrumModel = spectrumModel;
		this.aMap = this.spectrumModel.getAlignmentMapping();
		this.dataType = this.spectrumModel.getDataType();

		stateCount = dataType.getStateCount();
		operation = SpectrumOperation.NONE;
		
		addModel(this.spectrumModel);
		
		preprocessLikelihoodAlignmentMap();
		getLogLikelihood();
		storeEverything();
		
	}
	

	private void preprocessLikelihoodAlignmentMap() {
		makeDirty();
		
		liS = new LikelihoodScaler(LOG_C);
		
		srpCount = aMap.getSrpCount();
		spectrumCount = spectrumModel.getSpectrumCount();
		spectrumLength = spectrumModel.getSpectrumLength();
		
		logLikelihood = Double.NEGATIVE_INFINITY;
		storedLogLikelihood = Double.NEGATIVE_INFINITY;
		
		allLogLikelihood = new double[srpCount][spectrumCount][spectrumLength];
		storedAllLogLikelihood = new double[srpCount][spectrumCount][spectrumLength];

		spectrumLogLikelihood = new double[srpCount][spectrumCount];
		storedSpectrumLogLikelihood = new double[srpCount][spectrumCount];
		
		eachSrpLikelihood = new double[srpCount];
		storedEachSrpLikelihood = new double[srpCount];

//		logBinomialDesnity = new HashMap<Integer, double[]>();
//		scaledLogBinomialDesnity = new HashMap<Integer, double[]>();
//		
//		allDelta = new int[srpCount];

		
//		int maxDist=0;
//		for (int s = 0; s < srpCount; s++) {
//			String srp = aMap.getSrpFragment(s);
//
//			int srLength = srp.length();
////			char[] srCharArray = reads.toCharArray();
////			double plambda = srLength * ERROR_RATE;
////			double logPlambda = Math.log(plambda);
//
//			int srLength1 = srLength+1;
//			maxDist = Math.max(maxDist, srLength1);
//			double[] logBinomD = new double[srLength1];
//			double[] scaledBinomD = new double[srLength1];
//			for (int i = 0; i < logBinomD.length; i++) {
//
//				logBinomD[i] = ArithmeticUtils.binomialCoefficientLog(srLength, i)+i*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE;
//				scaledBinomD[i] = liS.scale(logBinomD[i]); 
//			}
////			System.out.println(Arrays.toString(logBinomD));
////			logBinomialDesnity.put(srLength, logBinomD);
////			scaledLogBinomialDesnity.put(srLength, scaledBinomD);
//		}
		
		
//		counter = new int[maxDist];
//		System.out.println(maxDist);

//		int dist = LikelihoodUtils.hamDist(reads, subH);
//		double logProb = -plambda + dist * logPlambda - allFactorialLog[dist];
//		double logProb = ArithmeticUtils.binomialCoefficientLog(srLength, dist)+dist*LOG_ERROR_RATE+(srLength-dist)*LOG_ONE_MINUS_ERROR_RATE;

	}


    
	@Override
	public double getLogLikelihood(){

        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
        

	}
	
	protected double calculateLogLikelihood() {
		SpectrumOperationRecord operationReocrd = spectrumModel.getSpectrumOperationRecord();
		operation = operationReocrd.getOperation();
//		System.out.println(operation);
		double logLikelihood = Double.NEGATIVE_INFINITY;

		switch (operation) {
			case NONE:
				logLikelihood = calculateSrpLikelihoodFull();
				break;
//			case SWAPMULTI:
//				logLikelihood = calculateSrpLikelihoodMultiBasesSwap();
//				break;

			case SINGLE_DELTA:
//				System.out.println("single");
//				logLikelihood = calculateSrpLikelihoodFull();				
				logLikelihood = calculateSrpLikelihoodSingle();
				break;
			case COLUMN_DELTA:
				logLikelihood = calculateSrpLikelihoodColumn();
//				logLikelihood = calculateSrpLikelihoodFull();
				break;
			case MULTI_DELTA:
//				logLikelihood = calculateSrpLikelihoodFull();
				logLikelihood = calculateSrpLikelihoodMulti();
				break;
				//TODO make multiCalculation
				//logLikelihood = calculateSrpLikelihoodMulti();
//			case SWAPCOLUMN:
//				logLikelihood = calculateSrpLikelihoodSwapColumn();
//				break;
			case RECOMBINATION:
				logLikelihood = calculateSrpLikelihoodRecombination();
				break;
//			case PASS:
//				logLikelihood = storedLogLikelihood;
//				break;
			default:
				throw new IllegalArgumentException("Unknown operation type: "+operation);
	
			}
//	    double logLikelihood = calculateShoreReadLikelihood4();
//	    double logLikelihood = calculateShoreReadLikelihoodBinomialModel2();
	    
//	    timeTrial();
//		storeState();
		return logLikelihood;
	}
	
	private double calculateSrpLikelihoodFull() {

//		System.out.println("calculateSrpLikelihoodFull");
		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
//		int spectrumIndex = record.getSpectrumIndex(); 
//		int siteIndex = record.getAllSiteIndexs()[0];
//		System.out.println("In fullCaluclation\t"+spectrumIndex +"\t"+ siteIndex +"\t"+ 
//				Arrays.toString(spectrumModel.getSpectrum(spectrumIndex)
//				.getFrequencies(siteIndex)));
		for (int i = 0; i < srpCount; i++) {

			String srp = aMap.getSrpFragment(i);
			String fullSrp = aMap.getSrpFull(i);
			int start = aMap.getSrpStart(i);
			int end = aMap.getSrpEnd(i);
			
			double logLikelihood;
			liS.reset();
			for (int j = 0; j < spectrumCount; j++) {

				Spectrum spectrum = spectrumModel.getSpectrum(j);
				spectrumLogLikelihood[i][j] = 0;
				for (int k = start; k < end; k++) {
					double[] frequencies = spectrum.getFrequenciesAt(k);
					char srpChar = fullSrp.charAt(k);
					int state = dataType.getState(srpChar);
					if(state<stateCount){
						double likelihood = frequencies[state] * NOT_ERROR_RATE
								+ (1 - frequencies[state]) * ERROR_RATE;
						logLikelihood = Math.log(likelihood);
					}
					else{
						logLikelihood = LOG_ERROR_RATE;
					}
//					if(spectrumIndex>0 && (allLogLikelihood[i][j][k] != logLikelihood)){
//						System.out.println("DIFFLI:"+i +" "+j+" "+" "+k+
//								" "+ logLikelihood + " "+ allLogLikelihood[i][j][k] +" "+
//								storedAllLogLikelihood[i][j][k]);
//					}
					allLogLikelihood[i][j][k] = logLikelihood;
					spectrumLogLikelihood[i][j] += allLogLikelihood[i][j][k];
					
				}
//				System.out.println("SL: "+spectrumLogLikelihood);
				liS.scaleLogProb(spectrumLogLikelihood[i][j]);
				
//				int dist = LikelihoodUtils.Dist(start, end, srp, haplotypeModel.getAlignedSequenceString(j));
//				allDists[i][j]=dist;
				
//				liS.addScaledLogProb(logPD[dist]);
//				liS.scaleLogProb(logPD[dist]);
			}	
			
			eachSrpLikelihood[i] = liS.getLogLikelihood();
			

		}
		double logLikelihood = StatUtils.sum(eachSrpLikelihood);
		
//for (int i = 0; i < srpCount; i++) {
//	System.out.println(Arrays.toString(allDists[i]));
//}
//System.out.println("==");

//		double logLikelihood2 = calculateSrpLikelihoodFull();
//		if(logLikelihood != logLikelihood2){
//			System.out.println(logLikelihood +"\t"+ logLikelihood2 +"\t"+ (logLikelihood-logLikelihood2));
//		}
//		
//		System.out.println(logLikelihood);
		return logLikelihood;
	}

	

	private double calculateSrpLikelihoodSingle() {

		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int j = record.getSpectrumIndex(); 
		int k = record.getAllSiteIndexs()[0];

		Spectrum spectrum = spectrumModel.getSpectrum(j);
		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);
		
		for (int i : mapToSrp) {
			String fullSrp = aMap.getSrpFull(i);
			updateLikelihoodAtIJK(i, j, k, spectrum, fullSrp);
			updateEachSrpAt(i);
		}

		double logLikelihood = StatUtils.sum(eachSrpLikelihood);
		return logLikelihood;
	}

	private double calculateSrpLikelihoodMulti() {

		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();

		int[] siteIndexs = record.getAllSiteIndexs();
		int j= record.getSpectrumIndex(); 
		Spectrum spectrum = spectrumModel.getSpectrum(j);
		
		Set<Integer> allSrpPos = new HashSet<Integer>();
		for (int i = 0; i < siteIndexs.length; i++) {
			ArrayList<Integer> mapToSrp = aMap.getMapToSrp(siteIndexs[i]);
			allSrpPos.addAll(mapToSrp);
		}

		for (Integer integer : allSrpPos) {
			
			int i = integer;
			String fullSrp = aMap.getSrpFull(i);
		
			for (int s = 0; s < siteIndexs.length; s++) {
				int k = siteIndexs[s];
				if(allLogLikelihood[i][j][k]!=0){
					updateLikelihoodAtIJK(i, j, k, spectrum, fullSrp);
				}
			}
			updateEachSrpAt(i);


		}
		double logLikelihood = StatUtils.sum(eachSrpLikelihood);

		return logLikelihood;


	}
	private double calculateSrpLikelihoodColumn() {

		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int k = record.getColumnIndex();
		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);

		for (int i : mapToSrp) {
			String fullSrp = aMap.getSrpFull(i);
			
			for (int j = 0; j < spectrumCount; j++) {
				Spectrum spectrum = spectrumModel.getSpectrum(j);
				updateLikelihoodAtIJK(i, j, k, spectrum, fullSrp);

			}
			updateEachSrpAt(i);

		}
		double logLikelihood = StatUtils.sum(eachSrpLikelihood);
		return logLikelihood;
	}

	private double calculateSrpLikelihoodRecombination() {
		
		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int[] twoSpectrums = record.getRecombinationSpectrumIndex();
		int[] twoPositions = record.getRecombinationPositionIndex();
		
		Set<Integer> allSrpPos = new HashSet<Integer>();
		for (int i = twoPositions[0]; i < twoPositions[1]; i++) {
			ArrayList<Integer> mapToSrp = aMap.getMapToSrp(i);
			allSrpPos.addAll(mapToSrp);
		}
		for (Integer integer : allSrpPos) {
			int i = integer;
			String fullSrp = aMap.getSrpFull(i);

			for (int r = 0; r < TWO; r++) {
				int j = twoSpectrums[r];
				Spectrum spectrum = spectrumModel.getSpectrum(j);

				for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
					if(allLogLikelihood[i][j][k]!=0){
						updateLikelihoodAtIJK(i,j,k, spectrum, fullSrp);
					}
				}
			}
			updateEachSrpAt(i);

		}
		double logLikelihood = StatUtils.sum(eachSrpLikelihood);
		return logLikelihood;

	}


	private void updateEachSrpAt(int i) {
		liS.reset();
		for (int s = 0; s < spectrumCount; s++) {
			liS.scaleLogProb(spectrumLogLikelihood[i][s]);
		}
		eachSrpLikelihood[i] = liS.getLogLikelihood();
	}
	private void updateLikelihoodAtIJK(int i, int j, int k, Spectrum spectrum, String fullSrp) {
		
		char srpChar = fullSrp.charAt(k);
		int state = dataType.getState(srpChar);
		double[] frequencies = spectrum.getFrequenciesAt(k);
		if(state<stateCount){
			double likelihood = frequencies[state] * NOT_ERROR_RATE
					+ (1 - frequencies[state]) * ERROR_RATE;
			logLikelihood = Math.log(likelihood);
		}
		else{
			logLikelihood = LOG_ERROR_RATE;
		}

		spectrumLogLikelihood[i][j] -= allLogLikelihood[i][j][k]; 
		allLogLikelihood[i][j][k] = logLikelihood;
		spectrumLogLikelihood[i][j] += allLogLikelihood[i][j][k];
		
	}

	private void updateLikelihoodAtIJKOld(int i, int j, int k, Spectrum spectrum, String fullSrp) {
		
		char srpChar = fullSrp.charAt(k);
		int state = dataType.getState(srpChar);
		double[] frequencies = spectrum.getFrequenciesAt(k);
		if(state<stateCount){
			double likelihood = frequencies[state] * NOT_ERROR_RATE
					+ (1 - frequencies[state]) * ERROR_RATE;
			logLikelihood = Math.log(likelihood);
		}
		else{
			logLikelihood = LOG_ERROR_RATE;
		}

		spectrumLogLikelihood[i][j] -= allLogLikelihood[i][j][k]; 
		allLogLikelihood[i][j][k] = logLikelihood;
		spectrumLogLikelihood[i][j] += allLogLikelihood[i][j][k];
		
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == spectrumModel) {
            // treeModel has changed so recalculate the intervals
//            eventsKnown = false;
        }
//        System.err.println("Call handleModelChangedEvent in SpectrumAlignmentModel");
        makeDirty();
		
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		System.err.println("Call handleVariableChangedEvent in SpectrumAlignmentModel");
	}

//	@Deprecated
//	public void restoreStatePublicTest(){
//		restoreState();
//	}

	private Set<Integer> allSrpPos = new HashSet<Integer>(200);
	@Override
	public void storeState() {
		
		System.arraycopy(eachSrpLikelihood, 0, storedEachSrpLikelihood, 0, eachSrpLikelihood.length);
		storedLogLikelihood = logLikelihood;
//		System.arraycopy(storedEachLikelihood, 0, eachLikelihood, 0, eachLikelihood.length);

//		for (int i = 0; i < allDists.length; i++) {
//		    System.arraycopy(allDists[i], 0, storedAllDists[i], 0, allDists[0].length);
//		}
		SpectrumOperationRecord spectrumOperationRecord = spectrumModel.getSpectrumOperationRecord();
		SpectrumOperation operation = spectrumOperationRecord.getOperation();
//		int spectrumIndex;
//		int siteIndex; = spectrumOperationRecord.getAllSiteIndexs()[0];
		ArrayList<Integer> mapToSrp;// = aMap.getMapToSrp(siteIndex);
//		int[] siteIndexs;
		int count=0;
		int count2=0;
		int count3=0;
		switch (operation) {
		case NONE:
			break;

		case COLUMN_DELTA:
			int columnIndex = spectrumOperationRecord.getColumnIndex();
			mapToSrp = aMap.getMapToSrp(columnIndex);
			for (int sp = 0; sp < spectrumCount; sp++) {
				for (int i : mapToSrp) {
					storedAllLogLikelihood[i][sp][columnIndex] = allLogLikelihood[i][sp][columnIndex];
					storedSpectrumLogLikelihood[i][sp] = spectrumLogLikelihood[i][sp];
				}
			}

			break;
		case SINGLE_DELTA:
		case MULTI_DELTA:

			int spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			int[] siteIndexs = spectrumOperationRecord.getAllSiteIndexs();
			
			count = 0;
			allSrpPos.clear();
			for (int s = 0; s < siteIndexs.length; s++) {
				mapToSrp = aMap.getMapToSrp(siteIndexs[s]);
				allSrpPos.addAll(mapToSrp);
			}
			int j = spectrumIndex;
			for (int i : allSrpPos) {
//				for (int r = 0; r < TWO; r++) {
					for (int s = 0; s < siteIndexs.length; s++) {
						int k = siteIndexs[s];
//						int j = twoSpectrums[r];
						if(allLogLikelihood[i][j][k] != 0){
						storedAllLogLikelihood[i][j][k] = allLogLikelihood[i][j][k];
						storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
						count++;
						}
					}
//				}
			}
			
			count2 = 0;
			for (int s = 0; s < siteIndexs.length; s++) {
				mapToSrp = aMap.getMapToSrp(siteIndexs[s]);
				for (int i : mapToSrp) {
					storedAllLogLikelihood[i][spectrumIndex][siteIndexs[s]] = allLogLikelihood[i][spectrumIndex][siteIndexs[s]];
					storedSpectrumLogLikelihood[i][spectrumIndex] = spectrumLogLikelihood[i][spectrumIndex];
					count2++;
				}
			}
			System.out.println(count +"\t"+ count2 +"\t"+ (count-count2));
			
			break;
		case RECOMBINATION:
			int[] twoSpectrums = spectrumOperationRecord.getRecombinationSpectrumIndex();
			int[] twoPositions = spectrumOperationRecord.getRecombinationPositionIndex();
	
				count = 0;
	//			Set<Integer> allSrpPos = new HashSet<Integer>();
			long time1 = System.currentTimeMillis();
			for (int t = 0; t < 1e3; t++) {
				
				count=0;
				allSrpPos.clear();
				for (int i = twoPositions[0]; i < twoPositions[1]; i++) {
					mapToSrp = aMap.getMapToSrp(i);
					allSrpPos.addAll(mapToSrp);
				}
				for (int i : allSrpPos) {
					for (int r = 0; r < TWO; r++) {
						for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
							j = twoSpectrums[r];
							if(allLogLikelihood[i][j][k] != 0){
							storedAllLogLikelihood[i][j][k] = allLogLikelihood[i][j][k];
							storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
							count++;
							}
						}
					}
				}
			}
			long time2 = System.currentTimeMillis();
			System.out.print((time2-time1)+"\t" );
			
			time1 = System.currentTimeMillis();
			for (int t = 0; t < 1e3; t++) {
				
				count2 = 0;
				for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
					mapToSrp = aMap.getMapToSrp(k);
					for (int i : mapToSrp) {
						for (int r = 0; r < TWO; r++) {
							j = twoSpectrums[r];
							storedAllLogLikelihood[i][j][k] = allLogLikelihood[i][j][k];
							storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
							count2++;
						}
					}
				}

			}
			time2 = System.currentTimeMillis();

			System.out.print((time2-time1)+"\t" );
			time1 = System.currentTimeMillis();
			for (int t = 0; t < 1e3; t++) {
				
				int srpCount = aMap.getSrpCount();
				boolean[] srpSwitch = new boolean[srpCount];
				for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
					mapToSrp = aMap.getMapToSrp(k);
					for (int i : mapToSrp) {
						srpSwitch[i] = true;
					}
				}
				
				
				count3 = 0;
//				for (int i : allSrpPos) {
				for (int i = 0; i < srpSwitch.length; i++) {
					
					if(srpSwitch[i]){
					for (int r = 0; r < TWO; r++) {
						for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
							j = twoSpectrums[r];
							if(allLogLikelihood[i][j][k] != 0){
							storedAllLogLikelihood[i][j][k] = allLogLikelihood[i][j][k];
							storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
							count3++;
							}
						}
					}
					}
				}
				
//				for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
//					mapToSrp = aMap.getMapToSrp(k);
//					for (int i : mapToSrp) {
//						for (int r = 0; r < TWO; r++) {
//							j = twoSpectrums[r];
//							storedAllLogLikelihood[i][j][k] = allLogLikelihood[i][j][k];
//							storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
//							count3++;
//						}
//					}
//				}

			}
			time2 = System.currentTimeMillis();
			System.out.print((time2-time1)+"\t" );
			
			System.out.println(count +"\t"+ count2 +"\t"+ count3 +"\t"+ (count-count2));
			
			
			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: "+operation +"\tin"+ShortReadsSpectrumLikelihood.class.getSimpleName() );
//			storeEverything();
//			break;
		}


		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
//		System.out.println(spectrumIndex +"\t"+ siteIndex);
//		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(siteIndex);
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

	@Override
	protected void restoreState() {
		
//		System.err.println("SR likelihood restore: " + logLikelihood +"\t"+ storedLogLikelihood);
		logLikelihood = storedLogLikelihood;
		System.arraycopy(storedEachSrpLikelihood, 0, eachSrpLikelihood, 0, eachSrpLikelihood.length);
		
		SpectrumOperationRecord spectrumOperationRecord = spectrumModel.getSpectrumOperationRecord();
		SpectrumOperation operation = spectrumOperationRecord.getOperation();
//		int spectrumIndex;
//		int siteIndex = spectrumOperationRecord.getAllSiteIndexs()[0];
		ArrayList<Integer> mapToSrp;// = aMap.getMapToSrp(siteIndex);
//		int[] siteIndexs;
		
		switch (operation) {
		case NONE:
			break;
		case COLUMN_DELTA:
			int columnIndex = spectrumOperationRecord.getColumnIndex();
			mapToSrp = aMap.getMapToSrp(columnIndex);
			
			for (int sp = 0; sp < spectrumCount; sp++) {
				for (int i : mapToSrp) {
					allLogLikelihood[i][sp][columnIndex] = storedAllLogLikelihood[i][sp][columnIndex];
					spectrumLogLikelihood[i][sp] = storedSpectrumLogLikelihood[i][sp];
				}
			}

			break;
		case SINGLE_DELTA:
		case MULTI_DELTA:
			int spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			int[] siteIndexs = spectrumOperationRecord.getAllSiteIndexs();
			for (int s = 0; s < siteIndexs.length; s++) {
				mapToSrp = aMap.getMapToSrp(siteIndexs[s]);
				for (int i : mapToSrp) {
					allLogLikelihood[i][spectrumIndex][siteIndexs[s]] = storedAllLogLikelihood[i][spectrumIndex][siteIndexs[s]];
					spectrumLogLikelihood[i][spectrumIndex] = storedSpectrumLogLikelihood[i][spectrumIndex];
				}
			}
			break;
		case RECOMBINATION:
			int[] twoSpectrums = spectrumOperationRecord.getRecombinationSpectrumIndex();
			int[] twoPositions = spectrumOperationRecord.getRecombinationPositionIndex();
			
//			Set<Integer> allSrpPos = new HashSet<Integer>();
//			for (int i = twoPositions[0]; i < twoPositions[1]; i++) {
//				mapToSrp = aMap.getMapToSrp(i);
//				allSrpPos.addAll(mapToSrp);
//			}
//			
//			for (int i : allSrpPos) {
//				for (int r = 0; r < TWO; r++) {
//					for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
//						int j = twoSpectrums[r];
//						allLogLikelihood[i][j][k] = storedAllLogLikelihood[i][j][k];
//						spectrumLogLikelihood[i][j] = storedSpectrumLogLikelihood[i][j];
//					}
//				}
//			}

			
			for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
				mapToSrp = aMap.getMapToSrp(k);
				for (int i : mapToSrp) {
					for (int r = 0; r < TWO; r++) {
						int j = twoSpectrums[r];
						allLogLikelihood[i][j][k] = storedAllLogLikelihood[i][j][k];
						spectrumLogLikelihood[i][j] = storedSpectrumLogLikelihood[i][j];
					}
				}
			}
			
			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: "+operation +"\tin"+ShortReadsSpectrumLikelihood.class.getSimpleName() );

		}

//
//		
//		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
//		int spectrumIndex = record.getSpectrumIndex();
//		int siteIndex = record.getSiteIndex();
//
//		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(siteIndex);
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

	@Override
	protected void acceptState() {
		//Do nothing
	}
	private void restoreEverything(){
		
		for (int i = 0; i < srpCount; i++) {
//			System.err.println(i +"\t"+ allLogLikelihood[i][spectrumIndex][siteIndex] +"\t"+ storedAllLogLikelihood[i][spectrumIndex][siteIndex]);
//			storedAllLogLikelihood[i][spectrumIndex][siteIndex] = allLogLikelihood[i][spectrumIndex][siteIndex];
			System.arraycopy(storedSpectrumLogLikelihood[i],0, spectrumLogLikelihood[i], 0, spectrumCount);
			for (int j = 0; j < spectrumCount; j++) {
//				for (int k = 0; k < spectrumLength; k++) {
//					if(allLogLikelihood[i][j][k] != storedAllLogLikelihood[i][j][k]){
//						System.out.println("DIFFLI:"+i +" "+j+" "+" "+k+
//								" "+ allLogLikelihood[i][j][k] +" "+
//								storedAllLogLikelihood[i][j][k]);
//					}
//				}
				System.arraycopy(storedAllLogLikelihood[i][j], 0, allLogLikelihood[i][j], 0, spectrumLength);
			}
			
		}
	}
	private void storeEverything() {

		System.arraycopy(eachSrpLikelihood, 0, storedEachSrpLikelihood, 0, eachSrpLikelihood.length);
		storedLogLikelihood = logLikelihood;
//		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
//		int spectrumIndex = record.getSpectrumIndex();
//		int siteIndex = record.getSiteIndex();
//		System.out.println(spectrumIndex +"\t"+ siteIndex);
		for (int i = 0; i < srpCount; i++) {
//			System.err.println(i +"\t"+ allLogLikelihood[i][spectrumIndex][siteIndex] +"\t"+ storedAllLogLikelihood[i][spectrumIndex][siteIndex]);
//			storedAllLogLikelihood[i][spectrumIndex][siteIndex] = allLogLikelihood[i][spectrumIndex][siteIndex];
			System.arraycopy(spectrumLogLikelihood[i],0, storedSpectrumLogLikelihood[i], 0, spectrumCount);
			for (int j = 0; j < spectrumCount; j++) {
//				for (int k = 0; k < spectrumLength; k++) {
//					if(allLogLikelihood[i][j][k] != storedAllLogLikelihood[i][j][k]){
//						System.out.println("DIFFLI:"+i +" "+j+" "+" "+k+
//								" "+ allLogLikelihood[i][j][k] +" "+
//								storedAllLogLikelihood[i][j][k]);
//					}
//				}
				System.arraycopy(allLogLikelihood[i][j], 0, storedAllLogLikelihood[i][j], 0, spectrumLength);
			}
			
		}
	}

	public SpectrumOperation getOperation(){
		return operation;
	}
	
	@Override
	public Model getModel() {
		return this;
		
	}


	@Override
	public void makeDirty() {
        likelihoodKnown = false;
		
	}

	public double[] getEachLikelihood() {
		double[] copyOfValues = new double[eachSrpLikelihood.length];
        System.arraycopy(eachSrpLikelihood, 0, copyOfValues, 0, copyOfValues.length);
		return copyOfValues;
	}

	@Deprecated
		private double calculateSrpLikelihoodMulti2() {
	
			SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
			int spectrumIndex = record.getSpectrumIndex(); 
			int[] siteIndexs = record.getAllSiteIndexs();
	
			int j = spectrumIndex;
			Spectrum spectrum = spectrumModel.getSpectrum(j);
			Set<Integer> allSrp = new HashSet<Integer>();
			for (int s = 0; s < siteIndexs.length; s++) {
				int k = siteIndexs[s];
			
				ArrayList<Integer> mapToSrp = aMap.getMapToSrp(siteIndexs[s]);
				for (int i : mapToSrp) {
					allSrp.add(i);
					String fullSrp = aMap.getSrpFull(i);
					double logLikelihood;
	
	//					if(allLogLikelihood[i][j][k]!=0){
							double[] frequencies = spectrum.getFrequenciesAt(k);
							char srpChar = fullSrp.charAt(k);
							int state = dataType.getState(srpChar);
							if(state<stateCount){
								double likelihood = frequencies[state] * NOT_ERROR_RATE
										+ (1 - frequencies[state]) * ERROR_RATE;
								logLikelihood = Math.log(likelihood);
							}
							else{
								logLikelihood = LOG_ERROR_RATE;
							}
	
							spectrumLogLikelihood[i][j] -= allLogLikelihood[i][j][k]; 
							allLogLikelihood[i][j][k] = logLikelihood;
							spectrumLogLikelihood[i][j] += allLogLikelihood[i][j][k];
	//					}
	
					}
			}
			for (int i : allSrp) {
				
			
				liS.reset();
				for (j = 0; j < spectrumCount; j++) {
	//				double spectrumLogLikelihood = 0;
	//				for (k = start; k < end; k++) {
	//					spectrumLogLikelihood += allLogLikelihood[i][j][k];
	//				}
	//				spectrumLogLikelihood =	StatUtils.sum(allLogLikelihood[i][j], start, end);
					liS.scaleLogProb(spectrumLogLikelihood[i][j]);
				}
				eachSrpLikelihood[i] = liS.getLogLikelihood();
			}
	
			
			double logLikelihood = StatUtils.sum(eachSrpLikelihood);
	//		System.out.println(logLikelihood);
			return logLikelihood;
	
	
		}

}