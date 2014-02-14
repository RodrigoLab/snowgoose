package srp.spectrum.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.ArithmeticUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.google.common.primitives.Doubles;

import srp.dr.evolution.datatype.ShortReads;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.ShortRead;
import srp.haplotypes.likelihood.LikelihoodScaler;
import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import srp.spectrum.SpectrumOperationRecord;
import dr.app.beauti.util.NumberUtil;
import dr.evolution.datatype.DataType;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;
import dr.math.distributions.BetaDistribution;
import dr.util.Assert;

/*

Three methods to go through multiple sites
Method 1: Go through them, might be faster for small srp size/small number of bases change

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


Method 2:
get all unique sites with HashSet

	Set<Integer> allSrpPos = new HashSet<Integer>();
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
				}
			}
		}
	}
Method 3: boolean index array
				
	int srpCount = aMap.getSrpCount();
	boolean[] srpSwitch = new boolean[srpCount];
	Arrays.fill(srpSwitch, false);
	for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
		mapToSrp = aMap.getMapToSrp(k);
		for (int i : mapToSrp) {
			srpSwitch[i] = true;
		}
	}
				
	for (int i = 0; i < srpSwitch.length; i++) {
		if(srpSwitch[i]){
			for (int r = 0; r < TWO; r++) {
				j = twoSpectrums[r];
				for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
					if(allLogLikelihood[i][j][k] != 0){
						storedAllLogLikelihood[i][j][k] = allLogLikelihood[i][j][k];
						storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
					}
				}
			}
		}
	}
				

 */


public class ShortReadsSpectrumLikelihood  extends AbstractModelLikelihood {

	/**
	 * 
	 */
	private static final long serialVersionUID = 7438385718398999755L;

	private static final boolean DEBUG = false;
	
    public static final String SHORT_READ_LIKELIHOOD = "ShortReadSpectrumLikelihood";
	public static final String NAME = SHORT_READ_LIKELIHOOD;
	public static final double ERROR_RATE = 0.0107;
	public static final double NOT_ERROR_RATE = 1-ERROR_RATE;
	public static final double LOG_ERROR_RATE = Math.log(ERROR_RATE);
	public static final double LOG_NOT_ERROR_RATE = Math.log(NOT_ERROR_RATE);
//	public static final double LOG_ONE_MINUS_ERROR_RATE = Math.log(1-ERROR_RATE);
	public static final double C = 1e-200;
	public static final double LOG_C = Math.log(C);
	
	private static final double EVALUATION_TEST_THRESHOLD = 1e-8;
	private static final int TWO = 2;

	private static final int GAP_STATE = 17;  
	private final int AMBIGUOUS_STATE_COUNT;
	private final int STATE_COUNT;
	
	protected boolean likelihoodKnown;
	
	private int spectrumLength;
	private int spectrumCount;
	private int srpCount;

	private double logLikelihood;
	private double storedLogLikelihood;

	private double[] eachSrpLogLikelihood;
	private double[] storedEachSrpLogLikelihood;

	private double[] sumScaledSrpLogLikelihood;
	private double[] storedSumSrpLogLikelihood;

	
	private double[][] spectrumLogLikelihood;
	private double[][] storedSpectrumLogLikelihood;
	
	private double[][] scaledSpectrumLogLikelihood;
	private double[][] storedScaledSpectrumLogLikelihood;
//	private double[][] spectrumScaledLogLikelihood;
//	private double[][] storedSpectrumScaledLogLikelihood;
	
	@Deprecated
	private double[][][] allLogLikelihood;
	@Deprecated
	private double[][][] storedAllLogLikelihood;


	private AlignmentMapping aMap;
	private LikelihoodScaler liS;
	
	private SpectrumAlignmentModel spectrumModel;
	private DataType dataType;
	
	
	private boolean[] srpSwitch;
	private Set<Integer> allSrpPos;
	static double[] temp = new double[18];

	private MultiType type;
	
	double[][] allStateLogLikelihood = new double[spectrumLength][20]; 
	double[][] allStoredStateLogLikelihood = new double[spectrumLength][20];

	private int[][] matchCount;
	private int[][] storedMatchCount;
	

	private HashMap<Integer, double[]> logBinomialDesnity;
	private HashMap<Integer, double[]> scaledLogBinomialDesnity;
	private int[] srLength;
	

	public ShortReadsSpectrumLikelihood(SpectrumAlignmentModel spectrumModel){
		super(SHORT_READ_LIKELIHOOD);
		Arrays.fill(temp, 1e-10);


		type = MultiType.Array;
//		type = MultiType.Hash;
//		type = MultiType.All;
		
		////
		
		
		
		this.spectrumModel = spectrumModel;
		this.aMap = this.spectrumModel.getAlignmentMapping();
		
		this.dataType = this.spectrumModel.getDataType();
		STATE_COUNT = dataType.getStateCount();//4
		AMBIGUOUS_STATE_COUNT = dataType.getAmbiguousStateCount();//18

		likelihoodKnown = false;
		
		
		addModel(this.spectrumModel);
		
		preprocessLikelihoodAlignmentMap();
		getLogLikelihood();
		storeEverything();
		
	}
	

	private void preprocessLikelihoodAlignmentMap() {
//		makeDirty();
		
		liS = new LikelihoodScaler(LOG_C);
		
		srpCount = aMap.getSrpCount();
		spectrumCount = spectrumModel.getSpectrumCount();
		spectrumLength = spectrumModel.getSpectrumLength();
		
		this.srpSwitch = new boolean[srpCount];
		this.allSrpPos = new HashSet<Integer>();
		
		logLikelihood = Double.NEGATIVE_INFINITY;
		storedLogLikelihood = Double.NEGATIVE_INFINITY;
		
		allLogLikelihood = new double[srpCount][spectrumCount][spectrumLength];
		storedAllLogLikelihood = new double[srpCount][spectrumCount][spectrumLength];

		spectrumLogLikelihood = new double[srpCount][spectrumCount];
		storedSpectrumLogLikelihood = new double[srpCount][spectrumCount];
		
		scaledSpectrumLogLikelihood = new double[srpCount][spectrumCount];
		storedScaledSpectrumLogLikelihood = new double[srpCount][spectrumCount];
		
		matchCount = new int[srpCount][spectrumCount];
		storedMatchCount = new int[srpCount][spectrumCount];
		
		sumScaledSrpLogLikelihood = new double[srpCount];
		storedSumSrpLogLikelihood = new double[srpCount];

		
		eachSrpLogLikelihood = new double[srpCount];
		storedEachSrpLogLikelihood = new double[srpCount];

		allStateLogLikelihood = new double[spectrumLength][AMBIGUOUS_STATE_COUNT]; 
		allStoredStateLogLikelihood = new double[spectrumLength][AMBIGUOUS_STATE_COUNT];
		for (int i = 0; i < spectrumLength; i++) {
			Arrays.fill(allStateLogLikelihood[i], LOG_ERROR_RATE);
			Arrays.fill(allStoredStateLogLikelihood[i], LOG_ERROR_RATE);

		}
		
		srLength = new int[srpCount];

		logBinomialDesnity = new HashMap<Integer, double[]>();
		scaledLogBinomialDesnity = new HashMap<Integer, double[]>();
		
		//Binom
		for (int s = 0; s < srpCount; s++) {
			String srp = aMap.getSrpFragment(s);

			srLength[s] = srp.length();
//			char[] srCharArray = reads.toCharArray();
//			double plambda = srLength * ERROR_RATE;
//			double logPlambda = Math.log(plambda);

			int srLength1 = srLength[s]+1;

			double[] logBinomD = new double[srLength1];
			double[] scaledBinomD = new double[srLength1];
			for (int i = 0; i < logBinomD.length; i++) {

				logBinomD[i] = ArithmeticUtils.binomialCoefficientLog(srLength[s], i);
//				logBinomD[i] = i*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE;
				scaledBinomD[i] = liS.scale(logBinomD[i]); 
			}
//			System.out.println(Arrays.toString(logBinomD));
			logBinomialDesnity.put(srLength[s], logBinomD);
			scaledLogBinomialDesnity.put(srLength[s], scaledBinomD);
		}
		
		
	}


    
	@Override
	public double getLogLikelihood(){


//		long time1 = System.currentTimeMillis();
		
		
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        
//		long time2 = System.currentTimeMillis();
//		time += (time2-time1);
        
        return logLikelihood;


		
		

	}
	
	protected double calculateLogLikelihood() {
//		SpectrumOperationRecord operationReocrd = spectrumModel.getSpectrumOperationRecord();
		SpectrumOperation operation = spectrumModel.getSpectrumOperation();
		double logLikelihood = Double.NEGATIVE_INFINITY;
//System.err.println("calculateLikelihood\t"+operation);
//		operation = SpectrumOperation.FULL;
		switch (operation) {
			case NONE:
			case FULL:
				if(DEBUG){
					System.out.println("Calculate ShortReadLikelihood:\t"+operation);
				}
				logLikelihood = calculateSrpLikelihoodFull();
//				logLikelihood = calculateSrpLikelihoodFullMaster();
				break;
//			case SWAPMULTI:
//				logLikelihood = calculateSrpLikelihoodMultiBasesSwap();
//				break;

			case DELTA_SINGLE:
			case SWAP_SINGLE:
//				System.out.println("single");
//				logLikelihood = calculateSrpLikelihoodFull();				
				logLikelihood = calculateSrpLikelihoodSingle();
				break;
			case DELTA_COLUMN:
			case SWAP_COLUMN:
				logLikelihood = calculateSrpLikelihoodColumn();
//				logLikelihood = calculateSrpLikelihoodFull();
				break;
			case DELTA_MULTI:
			case SWAP_MULTI:
//				logLikelihood = calculateSrpLikelihoodFull();
				logLikelihood = calculateSrpLikelihoodMulti();
				
				break;
			case SWAP_SUBCOLUMN:
				logLikelihood = calculateSrpLikelihoodSwapSubColumn();
				break;
			case RECOMBINATION:
				logLikelihood = calculateSrpLikelihoodRecombination();
				break;

//			case MASTER:
//				logLikelihood = calculateSrpLikelihoodFullMaster()
//				break;
			default:
				throw new IllegalArgumentException("Unknown operation type: "+operation);
	
			}
//	    double logLikelihood = calculateShoreReadLikelihood4();
//	    double logLikelihood = calculateShoreReadLikelihoodBinomialModel2();
	    
//	    timeTrial();
//		storeState();
//		System.out.println("likelihood\t"+ logLikelihood);
//		if( (logLikelihood-this.logLikelihood)==0 ){
//		System.err.println("Delta: "+ (logLikelihood-this.logLikelihood) +"\t"+ logLikelihood +"\t"+ this.logLikelihood +"\t"+ operation);
////		System.out.println("Delta: "+ (logLikelihood-this.logLikelihood) +"\t"+ logLikelihood +"\t"+ this.logLikelihood +"\t"+ operation);
//		}
//		else{
//			System.out.println("Delta: "+ (logLikelihood-this.logLikelihood) +"\t"+ logLikelihood +"\t"+ this.logLikelihood +"\t"+ operation);
//		}
		return logLikelihood;
	}
	
	private double calculateSrpLikelihoodFull() {

//		System.out.println("calculateSrpLikelihoodFull");


		double[][][] allStateLogLikelihoodFull = new double[spectrumCount][spectrumLength][AMBIGUOUS_STATE_COUNT]; 
//		double[][][] allStateLogLikelihood2 = new double[spectrumCount][spectrumLength][AMBIGUOUS_STATE_COUNT];
		for (int j = 0; j < spectrumCount; j++) {
			Spectrum spectrum = spectrumModel.getSpectrum(j);
			for (int k = 0; k < spectrumLength; k++) {
//				allStateLogLikelihood2[j][k] = calculateStatesLogLikelihood(spectrum, k);
				SpectraParameter spectra = spectrum.getSpectra(k);
				calculateStatesLogLikelihood(spectra, allStateLogLikelihood[k]);
				System.arraycopy(allStateLogLikelihood[k], 0, allStateLogLikelihoodFull[j][k], 0, AMBIGUOUS_STATE_COUNT);
			}
		}
			
		
		double stateLogLikelihood;
		for (int i = 0; i < srpCount; i++) {

			String fullSrp = aMap.getSrpFull(i);
			int start = aMap.getSrpStart(i);
			int end = aMap.getSrpEnd(i);
			
			liS.reset();
			for (int j = 0; j < spectrumCount; j++) {

//				Spectrum spectrum = spectrumModel.getSpectrum(j);
				spectrumLogLikelihood[i][j] = 0;
				int match = 0;
				for (int k = start; k < end; k++) {
					int state = getStateAtK(fullSrp, k);
					stateLogLikelihood = allStateLogLikelihoodFull[j][k][state];
					
//					if(stateLogLikelihood != LOG_ERROR_RATE){
					if(stateLogLikelihood == LOG_NOT_ERROR_RATE){
						match++;
					}
					
					
//					double[] frequencies = spectrum.getFrequenciesAt(k);
//					if(state<STATE_COUNT){
//						double likelihood = frequencies[state] * NOT_ERROR_RATE
//								+ (1 - frequencies[state]) * ERROR_RATE;
//						stateLogLikelihood = Math.log(likelihood);
//					}
//					else{
//						stateLogLikelihood = LOG_ERROR_RATE;
//					}

					allLogLikelihood[i][j][k] = stateLogLikelihood;
					spectrumLogLikelihood[i][j] += stateLogLikelihood;
					
				}
				matchCount[i][j] = match;
				double temp = logBinomialDesnity.get(srLength[i])[match];
//				spectrumLogLikelihood[i][j] +=ArithmeticUtils.binomialCoefficientLog((end-start), match);
				spectrumLogLikelihood[i][j] += temp;
//				spectrumScaledLogLikelihood[i][j] = liS.scale(spectrumLogLikelihood[i][j]);
				scaledSpectrumLogLikelihood[i][j] = liS.scale(spectrumLogLikelihood[i][j]);
				liS.add(scaledSpectrumLogLikelihood[i][j]);
//				System.out.println(spectrumLogLikelihood[i][j] +"\t"+ scaledSpectrumLogLikelihood[i][j] +"\t"+ liS.getLogLikelihood());
			}	
			sumScaledSrpLogLikelihood[i] = liS.getSumScaledLikelihood();
			eachSrpLogLikelihood[i] = liS.getLogLikelihood();
//			System.out.println(eachSrpLogLikelihood[i]);
		}
		
//		double logLikelihood = liS.sumLogLikelihood(sumScaledSrpLogLikelihood);
		double totalLogLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return totalLogLikelihood;
	}


	private double calculateSrpLikelihoodSingleBinom() {


		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int j = record.getSpectrumIndex(); 
		int k = record.getAllSiteIndexs()[0];
		int[] siteIndexs = record.getAllSiteIndexs();
//
		SpectraParameter spectra = spectrumModel.getSpectrum(j).getSpectra(k);
		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);
		calculateStatesLogLikelihood(spectra, allStateLogLikelihood[k]);
		calculateStoredStatesLogLikelihood(spectra, allStoredStateLogLikelihood[k]);
//		double[] stateLogLikelihood = calculateStatesLogLikelihood(spectrum, k);
//		double[] storedAllStateLogLikelihood = calculateStoredStatesLogLikelihood(spectrum.getSpectra(k));
		double totalLogLikelihood = this.logLikelihood;
//		System.out.println(j +"\t"+ k);
//		System.out.println(totalLogLikelihood);
		for (int i : mapToSrp) {
			String fullSrp = aMap.getSrpFull(i);
			int state = getStateAtK(fullSrp, k);
			
			totalLogLikelihood = updateLikelihoodAtIJK(i, j, k, state,
					allStateLogLikelihood[k], allStoredStateLogLikelihood[k],
					totalLogLikelihood);
		}
//		for (int i : mapToSrp) {
//			String fullSrp = aMap.getSrpFull(i);
//			char srpChar = fullSrp.charAt(k);
//			int state = dataType.getState(srpChar);
//				totalLogLikelihood = updateLikelihoodAtIJK(i, j, k, state, 
//						allStateLogLikelihood[k], allStoredStateLogLikelihood[k], totalLogLikelihood);
//		}
		return totalLogLikelihood;
	}
	
	private static int calculateDeltaDist(int srpChar, int newChar, int oldChar){//, boolean isHapEqualNew){
		
		int deltaDist = 0;
	
		if(newChar!= oldChar){ // if(newChar!= oldChar && isHapEqualNew)
			if (srpChar==newChar){
				deltaDist = -1;
			}
			else if(srpChar==oldChar){
				deltaDist = 1;
			}
		}
		return deltaDist;
		
	}
	
	private double calculateSrpLikelihoodSingle() {


		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int j = record.getSpectrumIndex(); 
		int k = record.getAllSiteIndexs()[0];
		int[] siteIndexs = record.getAllSiteIndexs();
//
		SpectraParameter spectra = spectrumModel.getSpectrum(j).getSpectra(k);
		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);
		calculateStatesLogLikelihood(spectra, allStateLogLikelihood[k]);
		calculateStoredStatesLogLikelihood(spectra, allStoredStateLogLikelihood[k]);
//		double[] stateLogLikelihood = calculateStatesLogLikelihood(spectrum, k);
//		double[] storedAllStateLogLikelihood = calculateStoredStatesLogLikelihood(spectrum.getSpectra(k));
		double totalLogLikelihood = this.logLikelihood;
//		System.out.println(j +"\t"+ k);
//		System.out.println(totalLogLikelihood);
		for (int i : mapToSrp) {
			String fullSrp = aMap.getSrpFull(i);
			int state = getStateAtK(fullSrp, k);

			totalLogLikelihood = updateLikelihoodAtIJK(i, j, k, state,
					allStateLogLikelihood[k], allStoredStateLogLikelihood[k],
					totalLogLikelihood);
		}
//		for (int i : mapToSrp) {
//			String fullSrp = aMap.getSrpFull(i);
//			char srpChar = fullSrp.charAt(k);
//			int state = dataType.getState(srpChar);
//				totalLogLikelihood = updateLikelihoodAtIJK(i, j, k, state, 
//						allStateLogLikelihood[k], allStoredStateLogLikelihood[k], totalLogLikelihood);
//		}
		return totalLogLikelihood;
	}
	

	private double calculateSrpLikelihoodMulti() {
		
		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();

		int[] siteIndexs = record.getAllSiteIndexs();
		int j= record.getSpectrumIndex(); 
		Spectrum spectrum = spectrumModel.getSpectrum(j);
		double totalLogLikelihood = this.logLikelihood;
		
		
		for (int s = 0; s < siteIndexs.length; s++) {
			int k = siteIndexs[s];
			SpectraParameter spectra = spectrum.getSpectra(k);

			calculateStatesLogLikelihood(spectra, allStateLogLikelihood[s]);
			calculateStoredStatesLogLikelihood(spectra, allStoredStateLogLikelihood[s]);
			
		}
		if(type==MultiType.Array){

			for (int s : siteIndexs) {
				ArrayList<Integer> mapToSrp = aMap.getMapToSrp(s);
				for (int i : mapToSrp) {
					srpSwitch[i] = true;
				}
			}

			for (int i = 0; i < srpSwitch.length; i++) {
				if(srpSwitch[i]){

					totalLogLikelihood = updateLikelihoodAtIJ(i, j, siteIndexs,
							allStateLogLikelihood, allStoredStateLogLikelihood,
							totalLogLikelihood);

				}
			}
	
		}
		else if(type==MultiType.Hash){
			allSrpPos.clear();
			for (int s : siteIndexs) {
				ArrayList<Integer> mapToSrp = aMap.getMapToSrp(s);
				allSrpPos.addAll(mapToSrp);
			}
			for (int i : allSrpPos) {
					totalLogLikelihood = updateLikelihoodAtIJ(i, j, siteIndexs,
							allStateLogLikelihood, allStoredStateLogLikelihood,
							totalLogLikelihood);

//				}
			}
		}
		else if(type==MultiType.All){
			for (int s = 0; s < siteIndexs.length; s++) {
				int k = siteIndexs[s];
				ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);

				for (int i : mapToSrp) {
					String fullSrp = aMap.getSrpFull(i);
					int state = getStateAtK(fullSrp, k);
					totalLogLikelihood = updateLikelihoodAtIJK(i, j, k,
							state, allStateLogLikelihood[s],
							allStoredStateLogLikelihood[s],
							totalLogLikelihood);
				}
			}
			
		}
		
		

		return totalLogLikelihood;

	}


	private double calculateSrpLikelihoodColumn() {

		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int k = record.getColumnIndex();
		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);

		double totalLogLikelihood = this.logLikelihood;
		
		for (int j = 0; j < spectrumCount; j++) {
			SpectraParameter spectra = spectrumModel.getSpectrum(j).getSpectra(k);
			
			calculateStatesLogLikelihood(spectra, allStateLogLikelihood[j]);
			calculateStoredStatesLogLikelihood(spectra, allStoredStateLogLikelihood[j]);
		}
		
		for (int i : mapToSrp) {
			
			String fullSrp = aMap.getSrpFull(i);
			int state = getStateAtK(fullSrp, k);
			
			for (int j = 0; j < spectrumCount; j++) {
				totalLogLikelihood = updateLikelihoodAtIJK(i, j, k, state,
						allStateLogLikelihood[j], allStoredStateLogLikelihood[j],
						totalLogLikelihood);

			}
		}
		return totalLogLikelihood;
	}
	
	private double calculateSrpLikelihoodSwapSubColumn() {

		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int k = record.getColumnIndex();
		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);
		int[] allSpectrumIndexs = record.getAllSpectrumIndexs();
		
		double totalLogLikelihood = this.logLikelihood;

//		for (int j = 0; j < spectrumCount; j++) {
		for (int j : allSpectrumIndexs) {
			SpectraParameter spectra = spectrumModel.getSpectrum(j).getSpectra(k);
			
			calculateStatesLogLikelihood(spectra, allStateLogLikelihood[j]);
			calculateStoredStatesLogLikelihood(spectra, allStoredStateLogLikelihood[j]);
		}
		
		for (int i : mapToSrp) {
			String fullSrp = aMap.getSrpFull(i);
			int state = getStateAtK(fullSrp, k);
//			for (int j = 0; j < spectrumCount; j++) {
			for (int j : allSpectrumIndexs) {
				totalLogLikelihood = updateLikelihoodAtIJK(i, j, k, state,
						allStateLogLikelihood[j], allStoredStateLogLikelihood[j],
						totalLogLikelihood);

			}


		}
		return logLikelihood;
	}

	
	private double calculateSrpLikelihoodRecombination() {
		
		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int[] twoSpectrums = record.getRecombinationSpectrumIndex();
		int[] twoPositions = record.getRecombinationPositionIndex();

//		Spectrum[] spectrums = new Spectrum[] {
//				spectrumModel.getSpectrum(twoSpectrums[0]),
//				spectrumModel.getSpectrum(twoSpectrums[1]) };
	
		int j0 = twoSpectrums[0];
		int j1 = twoSpectrums[1];
		int length = twoPositions[1] - twoPositions[0];
		
		Spectrum spectrum0 = spectrumModel.getSpectrum(j0); 
		Spectrum spectrum1 = spectrumModel.getSpectrum(j1);
		
		int[] siteIndexs = new int[length];
		for (int k = twoPositions[0], s=0; k < twoPositions[1]; k++, s++) {
			siteIndexs[s] = k;
			ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				srpSwitch[i] = true;
			}
		}
		
		for (int s = 0; s < siteIndexs.length; s++) {
			
			int k = siteIndexs[s];
			SpectraParameter spectra0 = spectrum0.getSpectra(k);
			SpectraParameter spectra1 = spectrum1.getSpectra(k);
//			allStoredStateLogLikelihood0[s] = calculateStoredStatesLogLikelihood(spectra0);
//			allStoredStateLogLikelihood1[s] = calculateStoredStatesLogLikelihood(spectra1);
			calculateStoredStatesLogLikelihood(spectra0, allStateLogLikelihood[s]);
			calculateStoredStatesLogLikelihood(spectra1, allStoredStateLogLikelihood[s]);
			

//			allStateLogLikelihood0[s] = calculateStatesLogLikelihood(spectrum0, k);
//			allStateLogLikelihood1[s] = calculateStatesLogLikelihood(spectrum1, k);

		}
		
		
//		Set<Integer> allSrpPos = new HashSet<Integer>();
//		allSrpPos.clear();
//		for (int i = twoPositions[0]; i < twoPositions[1]; i++) {
//			ArrayList<Integer> mapToSrp = aMap.getMapToSrp(i);
//			allSrpPos.addAll(mapToSrp);
//		}

		double totalLikelihood = this.logLikelihood;
		for (int i = 0; i < srpSwitch.length; i++) {
			if(srpSwitch[i]){
	
//				String fullSrp = aMap.getSrpFull(i);

				totalLikelihood = updateLikelihoodAtIJ(i, j0, siteIndexs,
//				totalLikelihood = updateLikelihoodAtIJ(i, j0, twoPositions[0], length,
						allStoredStateLogLikelihood,
						allStateLogLikelihood,
						totalLikelihood);
				totalLikelihood = updateLikelihoodAtIJ(i, j1, siteIndexs,
////				totalLikelihood = updateLikelihoodAtIJ(i, j1, twoPositions[0], length ,
						allStateLogLikelihood,
						allStoredStateLogLikelihood,
						totalLikelihood);
//				for (int s = 0; s < length; s++) {
//					int k = twoPositions[0]+s;
//
//					char srpChar = fullSrp.charAt(k);
//					int state = dataType.getState(srpChar);
//					
//					if(state!=18){
//						totalLikelihood = updateLikelihoodAtIJK(i, j0, k,
//								state, allStoredStateLogLikelihood1[s],
//								allStoredStateLogLikelihood0[s],
//								totalLikelihood);
//						totalLikelihood = updateLikelihoodAtIJK(i, j1, k,
//								state, allStoredStateLogLikelihood0[s],
//								allStoredStateLogLikelihood1[s],
//								totalLikelihood);
//					}
//
//				}
			}

		}
	
//		totalLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return totalLikelihood;

	}

	@Deprecated
	private double calculateSrpLikelihoodRecombination_swap() {
		
		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int[] twoSpectrums = record.getRecombinationSpectrumIndex();
		int[] twoPositions = record.getRecombinationPositionIndex();

//		Spectrum[] spectrums = new Spectrum[] {
//				spectrumModel.getSpectrum(twoSpectrums[0]),
//				spectrumModel.getSpectrum(twoSpectrums[1]) };
	
		int j0 = twoSpectrums[0];
		int j1 = twoSpectrums[1];
		
		for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
			
			ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);
//			System.out.println("Site: "+k +"\t"+ mapToSrp.size());
			for (int i : mapToSrp) {
				srpSwitch[i] = true;
			}
		}
		
		Set<Integer> allSrpPos = new HashSet<Integer>();
		allSrpPos.clear();
		for (int i = twoPositions[0]; i < twoPositions[1]; i++) {
			ArrayList<Integer> mapToSrp = aMap.getMapToSrp(i);
			allSrpPos.addAll(mapToSrp);
		}


		for (int i = 0; i < srpSwitch.length; i++) {
			if(srpSwitch[i]){
	
//				String fullSrp = aMap.getSrpFull(i);
//				for (int r = 0; r < TWO; r++) {
				double LL0 = storedSpectrumLogLikelihood[i][j0];
				double LL1 = storedSpectrumLogLikelihood[i][j1];
				for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
					if(storedAllLogLikelihood[i][j0][k]!=0){

						double L0 = storedAllLogLikelihood[i][j0][k];
						double L1 = storedAllLogLikelihood[i][j1][k];

						LL0 += (-L0 + L1);
						LL1 += (-L1 + L0);

						allLogLikelihood[i][j1][k] = L0;
						allLogLikelihood[i][j0][k] = L1;
					}

				}
				spectrumLogLikelihood[i][j0] = LL0;
				spectrumLogLikelihood[i][j1] = LL1;
				updateEachSrpAt(i);
			}

		}
	
		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return logLikelihood;

	}

	private double[] calculateStoredStatesLogLikelihood(SpectraParameter spectra, 
			double[] statesLogLikelihood) {
		for (int state = 0; state < STATE_COUNT; state++) {
			double frequency = spectra.getStoredFrequency(state);
			statesLogLikelihood[state] = caluclateStateLogLikelihood(frequency);
		}
		return statesLogLikelihood;
	}

	private double[] calculateStatesLogLikelihood(SpectraParameter spectra, 
			double[] statesLogLikelihood) {
		for (int state = 0; state < STATE_COUNT; state++) {
			double frequency = spectra.getFrequency(state);
			statesLogLikelihood[state] = caluclateStateLogLikelihood(frequency);
		}
		
		return statesLogLikelihood;
	}



	private int getStateAtK(String fullSrp, int k) {
		char srpChar = fullSrp.charAt(k);
		int state = dataType.getState(srpChar);
	
		return state;
	}


	private final double caluclateStateLogLikelihood(double frequency) {
		double logLikelihood = Math.log(frequency * NOT_ERROR_RATE
				+ (1 - frequency) * ERROR_RATE);
		return logLikelihood;
	}
	
	//XXX: beta: don't think this will work
//	static BetaDistribution betaD = new BetaDistribution(296.79, 3.21);
//	static double high = betaD.logPdf(0.98);
//	static double low = LOG_ERROR_RATE;//betaD.logPdf(0.8);
//	private final double caluclateStateLogLikelihood(double frequency) {
//		double logLikelihood = Math.log(frequency * NOT_ERROR_RATE
//				+ (1 - frequency) * ERROR_RATE);
//		System.out.println(high +"\t"+ low +"\t"+ logLikelihood);
//		logLikelihood = high;
//		if(frequency!=1){
//			logLikelihood = low;
//		}
////		org.apache.commons.math3.distribution.BetaDistribution b
////		System.out.println(logLikelihood);
//		return logLikelihood;
//	}


	private double updateLikelihoodAtIJK(int i, int j, int k, int state,
				double[] statesLogLikelihood, double[] storedStatesLogLikelihood,
				double totalLogLikelihood) {
		
			double stateLn= statesLogLikelihood[state];
			double storedStateLn = storedStatesLogLikelihood[state];

			if(storedStateLn != stateLn){
				double[] ds = logBinomialDesnity.get(srLength[i]);
//				System.out.println(storedStateLn +"\t"+ stateLn +"\t"+ LOG_ERROR_RATE);
				totalLogLikelihood -= eachSrpLogLikelihood[i];
	
				liS.setsumScaledLikelihood(sumScaledSrpLogLikelihood[i]);
				liS.minus(scaledSpectrumLogLikelihood[i][j]);
	//			liS.minusScaleLogProb( spectrumLogLikelihood[i][j]);
	
				spectrumLogLikelihood[i][j] -= ds[matchCount[i][j]];
						if(stateLn == LOG_NOT_ERROR_RATE){
							matchCount[i][j]+=1;
						}
						else{
							matchCount[i][j]-=1;
						}						
				spectrumLogLikelihood[i][j] -= storedStateLn; 
				spectrumLogLikelihood[i][j] += stateLn;
				try {
					int a = matchCount[i][j];
					double b = ds[matchCount[i][j]];
				} catch (Exception e) {
					System.out.println(i +"\t"+ j +"\t"+ matchCount[i][j] +"\t"+ storedMatchCount[i][j] +"\t"+ ds.length);
					System.out.println(aMap.getSrpFull(i).charAt(k) +"\t"+ Arrays.toString(spectrumModel.getSpectrum(j).getSpectra(k).getFrequencies()));
				}
				spectrumLogLikelihood[i][j] += ds[matchCount[i][j]];
				
				scaledSpectrumLogLikelihood[i][j] = liS.scale(spectrumLogLikelihood[i][j]);
				liS.add(scaledSpectrumLogLikelihood[i][j]);
	//			liS.addScaleLogProb(spectrumLogLikelihood[i][j]);
				sumScaledSrpLogLikelihood[i] = liS.getSumScaledLikelihood();
				eachSrpLogLikelihood[i] = liS.getLogLikelihood();
				
				totalLogLikelihood += eachSrpLogLikelihood[i];
	
			}
			return totalLogLikelihood;
		}


	private double updateLikelihoodAtIJ(int i, int j, int[] siteIndexs, 
				double[][] allStateLogLikelihood, double[][] storedAllStateLogLikelihood, 
				double totalLogLikelihood) {

		String fullSrp = aMap.getSrpFull(i);
		
		totalLogLikelihood -= eachSrpLogLikelihood[i];
		liS.setsumScaledLikelihood(sumScaledSrpLogLikelihood[i]);
		liS.minus( scaledSpectrumLogLikelihood[i][j]);
		
		for (int s = 0; s < siteIndexs.length; s++) {
			int k = siteIndexs[s];
			int state = getStateAtK(fullSrp, k);

			double stateLn= allStateLogLikelihood[s][state];
			double storedStateLn = storedAllStateLogLikelihood[s][state];
			if(storedStateLn != stateLn){
				spectrumLogLikelihood[i][j] -= storedStateLn; 
				spectrumLogLikelihood[i][j] += stateLn;
			}

		}
		scaledSpectrumLogLikelihood[i][j] = liS.scale(spectrumLogLikelihood[i][j]);
		liS.add(scaledSpectrumLogLikelihood[i][j]);
//		liS.addScaleLogProb(spectrumLogLikelihood[i][j]);
		sumScaledSrpLogLikelihood[i] = liS.getSumScaledLikelihood();
		eachSrpLogLikelihood[i] = liS.getLogLikelihood();
		totalLogLikelihood += eachSrpLogLikelihood[i];

		return totalLogLikelihood;
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == spectrumModel) {
            // treeModel has changed so recalculate the intervals
//            eventsKnown = false;
        }
        else{
        	System.err.println("Call handleModelChangedEvent in ShortReadSpectrumLikelihood" +"\t"+ model.getModelName());
        }
        likelihoodKnown = false;
//        System.err.println("Call handleModelChangedEvent in ShortReadSpectrumLikelihood" +"\t"+ model.getModelName());
//        makeDirty();
        
		
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		System.err.println("Call handleVariableChangedEvent in SpectrumAlignmentModel");
	}

	@Override
	protected void storeState() {
//long time1 = System.currentTimeMillis();

//		System.arraycopy(eachSrpLogLikelihood, 0, storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		storedLogLikelihood = logLikelihood;
//		storeEverything();
		SpectrumOperationRecord spectrumOperationRecord = spectrumModel.getSpectrumOperationRecord();
		SpectrumOperation operation = spectrumOperationRecord.getOperation();
//		int spectrumIndex;
//		int siteIndex; = spectrumOperationRecord.getAllSiteIndexs()[0];
		ArrayList<Integer> mapToSrp;// = aMap.getMapToSrp(siteIndex);

		int j;
		int k;
		switch (operation) {
		case NONE:
			if(DEBUG){
				System.out.println("StoreState in ShortReadsSpectrumLikelihood:\t"+operation);
			}
			break;
		case FULL:
			if(DEBUG){
				System.out.println("StoreState in ShortReadsSpectrumLikelihood:\t"+operation);
			}
			storeEverything();
			break;

		case DELTA_COLUMN:
		case SWAP_COLUMN:
		case SWAP_SUBCOLUMN:
			k= spectrumOperationRecord.getColumnIndex();
			mapToSrp = aMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
				storedSumSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i];
				for (j = 0; j < spectrumCount; j++) {
					storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
					storedScaledSpectrumLogLikelihood[i][j] = scaledSpectrumLogLikelihood[i][j];

				}
			}

			break;
		case DELTA_SINGLE:
		case SWAP_SINGLE:
			j = spectrumOperationRecord.getSpectrumIndex();
			k = spectrumOperationRecord.getAllSiteIndexs()[0];
			mapToSrp = aMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
				storedSumSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i];
				storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
				storedScaledSpectrumLogLikelihood[i][j] = scaledSpectrumLogLikelihood[i][j];
				storedMatchCount[i][j] = matchCount[i][j];
			}
			
		case DELTA_MULTI:
		case SWAP_MULTI:
			j = spectrumOperationRecord.getSpectrumIndex();
			int[] siteIndexs = spectrumOperationRecord.getAllSiteIndexs();

			if(type==MultiType.Array){
				for (int i = 0; i < srpSwitch.length; i++) {
					if (srpSwitch[i]) {
						storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
						storedSumSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i];
						storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
						storedScaledSpectrumLogLikelihood[i][j] = scaledSpectrumLogLikelihood[i][j];
					}
				}
			}
			else if(type==MultiType.Hash){
				for (int i : allSrpPos) {
					storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
					storedSumSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i];
					storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
					storedScaledSpectrumLogLikelihood[i][j] = scaledSpectrumLogLikelihood[i][j];
				}
			}
			else if(type==MultiType.All){
				for (int kk : siteIndexs) {
//				for (int s = 0; s < siteIndexs.length; s++) {
//					k = siteIndexs[s];
					mapToSrp = aMap.getMapToSrp(kk);
					for (int i : mapToSrp) {
						storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
						storedSumSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i];
						storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
						storedScaledSpectrumLogLikelihood[i][j] = scaledSpectrumLogLikelihood[i][j];

					}
				}
			}
						
			break;
		case RECOMBINATION:
			int[] twoSpectrums = spectrumOperationRecord.getRecombinationSpectrumIndex();
			int[] twoPositions = spectrumOperationRecord.getRecombinationPositionIndex();

			for (int i = 0; i < srpSwitch.length; i++) {
				if (srpSwitch[i]) {
					storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
					storedSumSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i];
					for (int r = 0; r < TWO; r++) {
						j = twoSpectrums[r];
						storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
						storedScaledSpectrumLogLikelihood[i][j] = scaledSpectrumLogLikelihood[i][j];

//						for ( k = twoPositions[0]; k < twoPositions[1]; k++) {
//							if (allLogLikelihood[i][j][k] != 0) {
//								storedAllLogLikelihood[i][j][k] = allLogLikelihood[i][j][k];
//								
//							}
//						}
					}
					srpSwitch[i] = false;
				}
			}
			
			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: "+operation +"\tin"+ShortReadsSpectrumLikelihood.class.getSimpleName() );
//			storeEverything();
//			break;
			
		}


//		long time2 = System.currentTimeMillis();
//		time += (time2-time1);
		
		
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
//		long time1 = System.currentTimeMillis();
		
//		System.err.println("SR likelihood restore: " + logLikelihood +"\t"+ storedLogLikelihood);
		logLikelihood = storedLogLikelihood;
//		restoreEverything();
//		System.arraycopy(storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		
		SpectrumOperationRecord spectrumOperationRecord = spectrumModel.getSpectrumOperationRecord();
		SpectrumOperation operation = spectrumOperationRecord.getOperation();
//		int spectrumIndex;
//		int siteIndex = spectrumOperationRecord.getAllSiteIndexs()[0];
		ArrayList<Integer> mapToSrp;// = aMap.getMapToSrp(siteIndex);
//		int[] siteIndexs;
		int j;
		int k;
		switch (operation) {
		case NONE:
			if(DEBUG){
				System.out.println("RestoreState in ShortReadsSpectrumLikelihood:\t"+operation);
			}
			break;
		case FULL:
			if(DEBUG){
				System.out.println("RestoreState in ShortReadsSpectrumLikelihood:\t"+operation);
			}
			restoreEverything();
			break;
		case DELTA_COLUMN:
		case SWAP_COLUMN:
		case SWAP_SUBCOLUMN:
			k = spectrumOperationRecord.getColumnIndex();
			mapToSrp = aMap.getMapToSrp(k);
			
			for (j = 0; j < spectrumCount; j++) {
				for (int i : mapToSrp) {
					eachSrpLogLikelihood[i] = storedEachSrpLogLikelihood[i];
					sumScaledSrpLogLikelihood[i] = storedSumSrpLogLikelihood[i]; 
					spectrumLogLikelihood[i][j] = storedSpectrumLogLikelihood[i][j];
					scaledSpectrumLogLikelihood[i][j] = storedScaledSpectrumLogLikelihood[i][j];
				}
			}

			break;
		case DELTA_SINGLE:
		case SWAP_SINGLE:
			j = spectrumOperationRecord.getSpectrumIndex();
			k = spectrumOperationRecord.getAllSiteIndexs()[0];
			mapToSrp = aMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				eachSrpLogLikelihood[i] = storedEachSrpLogLikelihood[i];
				sumScaledSrpLogLikelihood[i] = storedSumSrpLogLikelihood[i]; 
				spectrumLogLikelihood[i][j] = storedSpectrumLogLikelihood[i][j];
				scaledSpectrumLogLikelihood[i][j] = storedScaledSpectrumLogLikelihood[i][j];
				matchCount[i][j] = storedMatchCount[i][j];
			}
			
			
		case DELTA_MULTI:
		case SWAP_MULTI:
			j = spectrumOperationRecord.getSpectrumIndex();
			int[] siteIndexs = spectrumOperationRecord.getAllSiteIndexs();
			
			if(type==MultiType.Array){

				for (int i = 0; i < srpSwitch.length; i++) {
					if (srpSwitch[i]) {
						eachSrpLogLikelihood[i] = storedEachSrpLogLikelihood[i];
						sumScaledSrpLogLikelihood[i] = storedSumSrpLogLikelihood[i]; 
						spectrumLogLikelihood[i][j] = storedSpectrumLogLikelihood[i][j];
						scaledSpectrumLogLikelihood[i][j] = storedScaledSpectrumLogLikelihood[i][j];
					}
				}
			}

			else if(type==MultiType.Hash){
				for (int i : allSrpPos) {
					eachSrpLogLikelihood[i] = storedEachSrpLogLikelihood[i];
					sumScaledSrpLogLikelihood[i] = storedSumSrpLogLikelihood[i]; 
					spectrumLogLikelihood[i][j] = storedSpectrumLogLikelihood[i][j];
					scaledSpectrumLogLikelihood[i][j] = storedScaledSpectrumLogLikelihood[i][j];
				}
				
			}
			else if(type==MultiType.All){
				for (int kk : siteIndexs) {
					mapToSrp = aMap.getMapToSrp(kk);
					for (int i : mapToSrp) {
						eachSrpLogLikelihood[i] = storedEachSrpLogLikelihood[i];
						sumScaledSrpLogLikelihood[i] = storedSumSrpLogLikelihood[i]; 
						spectrumLogLikelihood[i][j] = storedSpectrumLogLikelihood[i][j];
						scaledSpectrumLogLikelihood[i][j] = storedScaledSpectrumLogLikelihood[i][j];
					}
				}

			}
			
			break;
		case RECOMBINATION:
			int[] twoSpectrums = spectrumOperationRecord.getRecombinationSpectrumIndex();
			int[] twoPositions = spectrumOperationRecord.getRecombinationPositionIndex();

			for (int i = 0; i < srpSwitch.length; i++) {
				if (srpSwitch[i]) {
					eachSrpLogLikelihood[i] = storedEachSrpLogLikelihood[i]; 
					sumScaledSrpLogLikelihood[i] = storedSumSrpLogLikelihood[i];
					for (int r = 0; r < TWO; r++) {
						j = twoSpectrums[r];
						spectrumLogLikelihood[i][j] = storedSpectrumLogLikelihood[i][j];
						scaledSpectrumLogLikelihood[i][j] = storedScaledSpectrumLogLikelihood[i][j];
					}
				}
			}

			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: "+operation +"\tin"+ShortReadsSpectrumLikelihood.class.getSimpleName() );

		}
//		long time2 = System.currentTimeMillis();
//		time += (time2-time1);
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
		
		System.arraycopy(storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		System.arraycopy(storedSumSrpLogLikelihood, 0, sumScaledSrpLogLikelihood, 0, sumScaledSrpLogLikelihood.length);
		for (int i = 0; i < srpCount; i++) {
//			System.err.println(i +"\t"+ allLogLikelihood[i][spectrumIndex][siteIndex] +"\t"+ storedAllLogLikelihood[i][spectrumIndex][siteIndex]);
//			storedAllLogLikelihood[i][spectrumIndex][siteIndex] = allLogLikelihood[i][spectrumIndex][siteIndex];
			System.arraycopy(storedSpectrumLogLikelihood[i], 0, spectrumLogLikelihood[i], 0, spectrumCount);
			System.arraycopy(storedScaledSpectrumLogLikelihood[i], 0, scaledSpectrumLogLikelihood[i], 0, spectrumCount);
			System.arraycopy(storedMatchCount[i], 0, matchCount[i], 0, spectrumCount);
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

		System.arraycopy(eachSrpLogLikelihood, 0, storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		System.arraycopy(sumScaledSrpLogLikelihood, 0, storedSumSrpLogLikelihood, 0, sumScaledSrpLogLikelihood.length);
//		storedLogLikelihood = logLikelihood;
//		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
//		int spectrumIndex = record.getSpectrumIndex();
//		int siteIndex = record.getSiteIndex();
//		System.out.println(spectrumIndex +"\t"+ siteIndex);
		for (int i = 0; i < srpCount; i++) {
//			System.err.println(i +"\t"+ allLogLikelihood[i][spectrumIndex][siteIndex] +"\t"+ storedAllLogLikelihood[i][spectrumIndex][siteIndex]);
//			storedAllLogLikelihood[i][spectrumIndex][siteIndex] = allLogLikelihood[i][spectrumIndex][siteIndex];
			System.arraycopy(spectrumLogLikelihood[i],0, storedSpectrumLogLikelihood[i], 0, spectrumCount);
			System.arraycopy(scaledSpectrumLogLikelihood[i],0, storedScaledSpectrumLogLikelihood[i], 0, spectrumCount);
			System.arraycopy(matchCount[i], 0, storedMatchCount[i], 0, spectrumCount);
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
		return spectrumModel.getSpectrumOperation();
	}
	
	@Override
	public Model getModel() {
		return this;
		
	}


	@Override
	public void makeDirty() {
//		System.err.println("make dirty");
		spectrumModel.resetSpectrumOperation();
        likelihoodKnown = false;
		
	}

	public double[] unittestMethodGetEachLikelihood() {
		double[] copyOfValues = new double[eachSrpLogLikelihood.length];
        System.arraycopy(eachSrpLogLikelihood, 0, copyOfValues, 0, copyOfValues.length);
		return copyOfValues;
	}


	@Override
	public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }


	public double calculateSrpLikelihoodFullMaster() {

		
		double logLikelihood = 0;
		double spectrumLogLikelihood = 0;
		double stateLogLikelihood = 0;
		
		for (int i = 0; i < srpCount; i++) {

			String fullSrp = aMap.getSrpFull(i);
			int start = aMap.getSrpStart(i);
			int end = aMap.getSrpEnd(i);
			
			liS.reset();
			for (int j = 0; j < spectrumCount; j++) {

				Spectrum spectrum = spectrumModel.getSpectrum(j);
				spectrumLogLikelihood = 0;
				for (int k = start; k < end; k++) {
					double[] frequencies = spectrum.getFrequenciesAt(k);
					int state = getStateAtK(fullSrp, k);
					if(state<STATE_COUNT){
						stateLogLikelihood = caluclateStateLogLikelihood(frequencies[state]);
//						double likelihood = frequencies[state] * NOT_ERROR_RATE
//								+ (1 - frequencies[state]) * ERROR_RATE;
//						stateLogLikelihood = Math.log(likelihood);
					}
					else{
						stateLogLikelihood = LOG_ERROR_RATE;
					}
					
					spectrumLogLikelihood += stateLogLikelihood;
					
				}
				liS.addLogProb(spectrumLogLikelihood);
				

			}	
			logLikelihood += liS.getLogLikelihood();
		}

//		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return logLikelihood;
	}



	public double calculateSrpLikelihoodFullMasterBinom() {
		//TODO unittest!!

		
		double logLikelihood = 0;
		double spectrumLogLikelihood = 0;
		double stateLogLikelihood = 0;
		
		for (int i = 0; i < srpCount; i++) {

			String fullSrp = aMap.getSrpFull(i);
			int start = aMap.getSrpStart(i);
			int end = aMap.getSrpEnd(i);
			
			liS.reset();
			for (int j = 0; j < spectrumCount; j++) {

				Spectrum spectrum = spectrumModel.getSpectrum(j);
				spectrumLogLikelihood = 0;
				for (int k = start; k < end; k++) {
					double[] frequencies = spectrum.getFrequenciesAt(k);
					int state = getStateAtK(fullSrp, k);
					if(state<STATE_COUNT){
						stateLogLikelihood = caluclateStateLogLikelihood(frequencies[state]);
//						double likelihood = frequencies[state] * NOT_ERROR_RATE
//								+ (1 - frequencies[state]) * ERROR_RATE;
//						stateLogLikelihood = Math.log(likelihood);
					}
					else{
						stateLogLikelihood = LOG_ERROR_RATE;
					}
					
					spectrumLogLikelihood += stateLogLikelihood;
					
				}
				liS.addLogProb(spectrumLogLikelihood);
				

			}	
			logLikelihood += liS.getLogLikelihood();
		}

//		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return logLikelihood;
	}


@Deprecated
private double calculateSrpLikelihoodSingle5() {


	SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
	int j = record.getSpectrumIndex(); 
	int k = record.getAllSiteIndexs()[0];

	Spectrum spectrum = spectrumModel.getSpectrum(j);
	ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);
	
	for (int i : mapToSrp) {
		String fullSrp = aMap.getSrpFull(i);
		updateLikelihoodAtIJK_getFreq(i, j, k, spectrum, fullSrp);
		updateEachSrpAt(i);
	}

	double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);

	return logLikelihood;
}


@Deprecated
private double calculateSrpLikelihoodMulti5() {

		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();

		int[] siteIndexs = record.getAllSiteIndexs();
		int j= record.getSpectrumIndex(); 
		Spectrum spectrum = spectrumModel.getSpectrum(j);
		
//		Arrays.fill(srpSwitch, true);
		for (int s : siteIndexs) {
			ArrayList<Integer> mapToSrp = aMap.getMapToSrp(s);
			for (int i : mapToSrp) {
				srpSwitch[i] = true;
			}
		}

	
//	Set<Integer> allSrpPos = new HashSet<Integer>();
//	allSrpPos.clear();
//	for (int s : siteIndexs) {
//
//		ArrayList<Integer> mapToSrp = aMap.getMapToSrp(s);
//		allSrpPos.addAll(mapToSrp);
//	}
//	System.out.println("Total:\t"+allSrpPos.size());


		for (int i = 0; i < srpSwitch.length; i++) {
			
			if(srpSwitch[i]){
				String fullSrp = aMap.getSrpFull(i);
			
				for (int s = 0; s < siteIndexs.length; s++) {
					int k = siteIndexs[s];
					if(allLogLikelihood[i][j][k]!=0){
						updateLikelihoodAtIJK_getFreq(i, j, k, spectrum, fullSrp);
					}
				}
				updateEachSrpAt(i);
			}

		}
		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);

		return logLikelihood;

}


@Deprecated
private double calculateSrpLikelihoodColumn5() {

	SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
	int k = record.getColumnIndex();
	ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);

	
	
	for (int i : mapToSrp) {
		String fullSrp = aMap.getSrpFull(i);
		
		for (int j = 0; j < spectrumCount; j++) {
			Spectrum spectrum = spectrumModel.getSpectrum(j);
			updateLikelihoodAtIJK_getFreq(i, j, k, spectrum, fullSrp);

		}
		updateEachSrpAt(i);

	}
	double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
	return logLikelihood;
}


@Deprecated
	private double calculateSrpLikelihoodRecombination_full() {
		
		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int[] twoSpectrums = record.getRecombinationSpectrumIndex();
		int[] twoPositions = record.getRecombinationPositionIndex();

		Spectrum[] spectrums = new Spectrum[] {
				spectrumModel.getSpectrum(twoSpectrums[0]),
				spectrumModel.getSpectrum(twoSpectrums[1]) };

		for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
			ArrayList<Integer> mapToSrp = aMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				srpSwitch[i] = true;
			}
		}
		for (int i = 0; i < srpSwitch.length; i++) {
			if(srpSwitch[i]){
				
				String fullSrp = aMap.getSrpFull(i);
//				for (int r = 0; r < TWO; r++) {
//					int j=twoSpectrums[r];	
					for (int k = twoPositions[0]; k < twoPositions[1]; k++) {
						
						int j=twoSpectrums[0];	
						if(allLogLikelihood[i][j][k]!=0){
							
							updateLikelihoodAtIJK_getFreq(i,j,k, spectrums[0], fullSrp);
							j=twoSpectrums[1];	
							updateLikelihoodAtIJK_getFreq(i,j,k, spectrums[1], fullSrp);
						}
					}
//				}
				updateEachSrpAt(i);
			}

		}

		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return logLikelihood;

	}


@Deprecated
	private void updateEachSrpAt(int i) {
		double temp = eachSrpLogLikelihood[i];
		liS.reset();
		for (int j = 0; j < spectrumCount; j++) {
			liS.addLogProb(spectrumLogLikelihood[i][j]);
		}
//		System.out.print("eachsrp" +"\t"+ eachSrpLikelihood[i]);
		eachSrpLogLikelihood[i] = liS.getLogLikelihood();
//		if(eachSrpLogLikelihood[i] != temp){
//			System.out.println("diff eachSrp: "+i +"\t"+ temp +"\t"+ eachSrpLogLikelihood[i]);
//		}
//		System.out.println(eachSrpLikelihood[i]);
	}


@Deprecated
	private void updateLikelihoodAtIJK_getFreq(int i, int j, int k, Spectrum spectrum, String fullSrp) {
		
		char srpChar = fullSrp.charAt(k);
		int state = dataType.getState(srpChar);
		double frequency = 0;
		double logLikelihood = LOG_ERROR_RATE;
		if(state<STATE_COUNT){
			frequency = spectrum.getFrequency(k, state);
			double likelihood = frequency * NOT_ERROR_RATE
					+ (1 - frequency) * ERROR_RATE;
			logLikelihood = Math.log(likelihood);
			caluclateStateLogLikelihood(frequency);
		}
//		System.out.println(i +"\t"+ j +"\t"+ k +"\t"+ state +"\t"+ frequency +"\t"+ logLikelihood);
//		System.out.print(spectrumLogLikelihood[i][j] +"\t" );
//		if(allLogLikelihood[i][j][k] != logLikelihood){
//			System.out.println("should be equal\t"+ allLogLikelihood[i][j][k] +"\t"+ logLikelihood);
//		}
		
		spectrumLogLikelihood[i][j] -= allLogLikelihood[i][j][k]; 
		allLogLikelihood[i][j][k] = logLikelihood;
		spectrumLogLikelihood[i][j] += allLogLikelihood[i][j][k];
//		System.out.println(spectrumLogLikelihood[i][j] +"\t"+ logLikelihood);
	}

}

enum MultiType{
	Array,
	Hash,
	All,;

};