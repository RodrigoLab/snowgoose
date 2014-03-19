package srp.spectrum.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
//import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.ArithmeticUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.carrotsearch.hppc.BitSet;
import com.carrotsearch.hppc.BitSetIterator;
import com.google.common.primitives.Doubles;

import srp.dr.evolution.datatype.ShortReads;
import srp.haplotypes.likelihood.LikelihoodScaler;
import srp.shortreads.AlignmentMapping;
import srp.shortreads.ShortRead;
import srp.shortreads.ShortReadMapping;
import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import srp.spectrum.SpectrumOperationRecord;
import srp.spectrum.likelihood.stateLikelihood.BetaMeanStateLikelihood;
import srp.spectrum.likelihood.stateLikelihood.BetaModeStateLikelihood;
import srp.spectrum.likelihood.stateLikelihood.ChisqStateLikelihood;
import srp.spectrum.likelihood.stateLikelihood.GTestStateLikelihood;
import srp.spectrum.likelihood.stateLikelihood.ProbabilityStateLikelihood;
import srp.spectrum.likelihood.stateLikelihood.StateLikelihood;
import dr.app.beauti.util.NumberUtil;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;
import dr.math.distributions.BetaDistribution;
import dr.math.distributions.ChiSquareDistribution;
import dr.math.distributions.GammaDistribution;
import dr.util.Assert;


public class ShortReadsSpectrumLikelihood  extends AbstractModelLikelihood {

	private static final long serialVersionUID = 7438385718398999755L;

	private static final boolean DEBUG = false;
	
    public static final String SHORT_READ_LIKELIHOOD = "ShortReadSpectrumLikelihood";
	public static final String NAME = SHORT_READ_LIKELIHOOD;
	public static final double ERROR_RATE = 0.0107;
	public static final double NOT_ERROR_RATE = 1-ERROR_RATE;
	public static final double LOG_ERROR_RATE = Math.log(ERROR_RATE);
	public static final double LOG_NOT_ERROR_RATE = Math.log(NOT_ERROR_RATE);
//	public static final double LOG_ONE_MINUS_ERROR_RATE = Math.log(1-ERROR_RATE);
	public static final double C = 1e-250;
	public static final double LOG_C = Math.log(C);

	public static final DataType DATA_TYPE = ShortReads.INSTANCE;
	public static final int STATE_COUNT = DATA_TYPE.getStateCount();//4
	public static final int AMBIGUOUS_STATE_COUNT = DATA_TYPE.getAmbiguousStateCount();//18
	
	private static final double EVALUATION_TEST_THRESHOLD = 1e-8;
	private static final int TWO = 2;
	private static final int GAP_STATE = 17;  
	

	
	private final double MIN_LOG_LIKELIHOOD;

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
	
	
	private double[] spectrumLogLikelihood;
	private double[] storedSpectrumLogLikelihood;
	
	private double[] scaledSpectrumLogLikelihood;
	private double[] storedScaledSpectrumLogLikelihood;
	
	
	@Deprecated private double[][] allStateLogLikelihood2D;
	@Deprecated private double[][] allStoredStateLogLikelihood2D;
	@Deprecated private double[][] spectrumLogLikelihood2D;
	@Deprecated private double[][] storedSpectrumLogLikelihood2D;
	@Deprecated private double[][] scaledSpectrumLogLikelihood2D;
	@Deprecated private double[][] storedScaledSpectrumLogLikelihood2D;
	
	private double[] allStateLogLikelihood;
	private double[] allStoredStateLogLikelihood;
	
	private double[] allStateLogLikelihoodIJK;
	private double[] allStoredStateLogLikelihoodIJK;
	
	private StateLikelihood stateLikelihood;
	private ShortReadMapping srpMap;
	private SpectrumAlignmentModel spectrumModel;
	private LikelihoodScaler liS;

	
	private int[] srpIndex ;
	private boolean[] srpSwitch;
	private Set<Integer> allSrpPos;

	private MultiType multiType;
//	private DistType distType;

	private String[] srpArray;
	private int[][] allSrpState2D;

	private int[][] mapToSrpArray;

	private BitSet bitSet;

	private int srpIndexCount;
	
//	public long timeStart = 0;//REMOVE
//	public long m1Time = 0;//REMOVE
//	public long m2Time = 0;//REMOVE
//	public long liTime  = 0;//REMOVE
//
//	double[][] m1State;// = new double[STATE_COUNT][spectrumLength];//REMOVE
//	double[][] m1StoreState;// = new double[STATE_COUNT][spectrumLength];//REMOVE
//	
//	double[][] m2State;// = new double[spectrumLength][STATE_COUNT];//REMOVE
//	double[][] m2StoreState;// = new double[spectrumLength][STATE_COUNT];//REMOVE
//	
//	double[] m3State;// = new double[STATE_COUNT][spectrumLength];//REMOVE
//	double[] m3StoreState;// = new double[STATE_COUNT][spectrumLength];//REMOVE

	public ShortReadsSpectrumLikelihood(SpectrumAlignmentModel spectrumModel, ShortReadMapping srpMap){
		this(spectrumModel, srpMap, DistType.flat);
		
	}
	
	public ShortReadsSpectrumLikelihood(SpectrumAlignmentModel spectrumModel, ShortReadMapping srpMap, DistType distType){
		super(SHORT_READ_LIKELIHOOD);
		this.spectrumModel = spectrumModel;
		this.srpMap = srpMap;

//		multiType = MultiType.Array;
		multiType = MultiType.BitSet;
		
//		type = MultiType.Hash;
//		type = MultiType.All;
//		distTypeCode = "flat";//"betaMean"  "betaMode" "gTest"
		setDistType(distType);
		MIN_LOG_LIKELIHOOD = stateLikelihood.caluclateStateLogLikelihood(SpectraParameter.MIN_FREQ);
		
		

		likelihoodKnown = false;
		
		addModel(this.spectrumModel);
		
		preprocessLikelihoodAlignmentMap();
		getLogLikelihood();
		
		
		storeEverything();
		
		for (int i = 0; i < srpCount; i++) {
			for (int j = 0; j < spectrumCount; j++) {
				Spectrum spectrum = spectrumModel.getSpectrum(j);
				for (int s = 0; s < spectrumLength; s++) {
					spectrum.getSpectra(s).storeState();
				}
			}
		}
	}
	

	private void preprocessLikelihoodAlignmentMap() {
//		makeDirty();
		
		liS = new LikelihoodScaler(LOG_C);
		
		srpCount = srpMap.getSrpCount();
		spectrumCount = spectrumModel.getSpectrumCount();
		spectrumLength = spectrumModel.getSpectrumLength();
		
		this.srpSwitch = new boolean[srpCount];
		this.allSrpPos = new HashSet<Integer>();
		this.srpIndex = new int[srpCount];
		
		logLikelihood = Double.NEGATIVE_INFINITY;
		storedLogLikelihood = Double.NEGATIVE_INFINITY;
//		
//		allLogLikelihood = new double[srpCount][spectrumCount][spectrumLength];
//		storedAllLogLikelihood = new double[srpCount][spectrumCount][spectrumLength];

		spectrumLogLikelihood2D = new double[srpCount][spectrumCount];
		storedSpectrumLogLikelihood2D = new double[srpCount][spectrumCount];
		
		scaledSpectrumLogLikelihood2D = new double[srpCount][spectrumCount];
		storedScaledSpectrumLogLikelihood2D = new double[srpCount][spectrumCount];
		
		spectrumLogLikelihood = new double[srpCount*spectrumCount];
		storedSpectrumLogLikelihood = new double[srpCount*spectrumCount];

		scaledSpectrumLogLikelihood = new double[srpCount*spectrumCount];
		storedScaledSpectrumLogLikelihood = new double[srpCount*spectrumCount];
		
		sumScaledSrpLogLikelihood = new double[srpCount];
		storedSumSrpLogLikelihood = new double[srpCount];

		
		eachSrpLogLikelihood = new double[srpCount];
		storedEachSrpLogLikelihood = new double[srpCount];

		allStateLogLikelihood = new double[spectrumLength*STATE_COUNT];
		allStoredStateLogLikelihood = new double[spectrumLength*STATE_COUNT];
		
		allStateLogLikelihood2D = new double[spectrumLength][AMBIGUOUS_STATE_COUNT]; 
		allStoredStateLogLikelihood2D = new double[spectrumLength][AMBIGUOUS_STATE_COUNT];
		
//		double errorLogLikelihood = stateLikelihood.caluclateStateLogLikelihood(SpectraParameter.MIN_FREQ);
		for (int i = 0; i < spectrumLength; i++) {
			Arrays.fill(allStateLogLikelihood2D[i], MIN_LOG_LIKELIHOOD);
			Arrays.fill(allStoredStateLogLikelihood2D[i], MIN_LOG_LIKELIHOOD);

		}
		
		srpArray = srpMap.getSrpArray();
		mapToSrpArray = srpMap.getMapToSrpArray();
		allSrpState2D = new int[srpArray.length][spectrumLength];

		for (int i = 0; i < srpArray.length; i++) {
			String srp = srpArray[i];
			for (int j = 0; j < spectrumLength; j++) {
				allSrpState2D[i][j] = getStateAtK(srp, j);

			}
		}
		
		bitSet = new BitSet(srpCount);
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
		
//		SpectrumOperationRecord operationReocrd = spectrumModel.getSpectrumOperationRecord();
		SpectrumOperation operation = spectrumModel.getSpectrumOperation();
		double logLikelihood = Double.NEGATIVE_INFINITY;
//System.err.println("calculateLikelihood\t"+operation);
//System.err.println("calculateLikelihood\t"+distType);
//		operation = SpectrumOperation.FULL;
		if(DEBUG){
			System.out.println("Calculate ShortReadLikelihood:\t"+operation);
		}
		switch (operation) {
			case NONE:
			
//				if(DEBUG){
//					System.out.println("Calculate ShortReadLikelihood:\t"+operation);
//				}
				logLikelihood = calculateSrpLikelihoodFull();
//				logLikelihood = calculateSrpLikelihoodFullMaster();
				break;
			case FULL:
//				if(DEBUG){
//					System.out.println("Calculate ShortReadLikelihood:\t"+operation);
//				}
				logLikelihood = calculateSrpLikelihoodFullMaster();
				break;
//			case SWAPMULTI:
//				logLikelihood = calculateSrpLikelihoodMultiBasesSwap();
//				break;

			case DELTA_SINGLE:
			case SWAP_SINGLE:
//				
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
				logLikelihood = calculateSrpLikelihoodSubColumn();
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


		double[][][] allStateLogLikelihoodFull3D = new double[spectrumCount][spectrumLength][AMBIGUOUS_STATE_COUNT]; 
		double[][] allStateLogLikelihoodFull2D = new double[spectrumCount][spectrumLength*STATE_COUNT];
//		double[][][] allStateLogLikelihood2 = new double[spectrumCount][spectrumLength][AMBIGUOUS_STATE_COUNT];
		for (int j = 0; j < spectrumCount; j++) {
			Spectrum spectrum = spectrumModel.getSpectrum(j);
			for (int k = 0; k < spectrumLength; k++) {
//				allStateLogLikelihood2[j][k] = calculateStatesLogLikelihood(spectrum, k);
				SpectraParameter spectra = spectrum.getSpectra(k);
				stateLikelihood.calculateStatesLogLikelihood2D(spectra, allStateLogLikelihood2D[k]);
				System.arraycopy(allStateLogLikelihood2D[k], 0, allStateLogLikelihoodFull3D[j][k], 0, AMBIGUOUS_STATE_COUNT );
				
				int kOffset = k*STATE_COUNT;
				stateLikelihood.calculateStatesLogLikelihood(spectra, kOffset, allStateLogLikelihood);
				
			}
			System.arraycopy(allStateLogLikelihood,0  , allStateLogLikelihoodFull2D[j], 0, allStateLogLikelihood.length );
		}
			
//		for (int s = 0; s < siteIndexs.length; s++) {
//			int k = siteIndexs[s];
//			SpectraParameter spectra = spectrum.getSpectra(k);
//			int kOffset = k*STATE_COUNT;
//			stateLikelihood.calculateStatesLogLikelihood(spectra, kOffset, allStateLogLikelihood);
//			stateLikelihood.calculateStoredStatesLogLikelihood(spectra, kOffset, allStoredStateLogLikelihood);
//			
//		}

		double stateLogLikelihood;
		for (int i = 0; i < srpCount; i++) {

//			String fullSrp = srpMap.getSrpFull(i);
			int start = srpMap.getSrpStart(i);
			int end = srpMap.getSrpEnd(i);
			
			liS.reset();
			for (int j = 0; j < spectrumCount; j++) {
				int offset = i*spectrumCount+j;
//				Spectrum spectrum = spectrumModel.getSpectrum(j);
				spectrumLogLikelihood2D[i][j] = 0;
				spectrumLogLikelihood[offset] = 0;
				for (int k = start; k < end; k++) {
//					int state = getStateAtK(fullSrp, k);
					int state = allSrpState2D[i][k];
					if(state<STATE_COUNT){
						stateLogLikelihood = allStateLogLikelihoodFull3D[j][k][state];
					}
					else{
						stateLogLikelihood = MIN_LOG_LIKELIHOOD;
					}
	
//					allLogLikelihood[i][j][k] = stateLogLikelihood;
					spectrumLogLikelihood2D[i][j] += stateLogLikelihood;
					spectrumLogLikelihood[offset] += stateLogLikelihood;
					
				}
				
				scaledSpectrumLogLikelihood[offset] = liS.scale(spectrumLogLikelihood[offset]);
				scaledSpectrumLogLikelihood2D[i][j] = liS.scale(spectrumLogLikelihood2D[i][j]);
//				liS.add(scaledSpectrumLogLikelihood2D[i][j]);
				liS.add(scaledSpectrumLogLikelihood[offset]);
//				if(i == 80){
//				System.out.println(i +"\t"+ j +"\t"+ spectrumLogLikelihood[i][j] +"\t"+ scaledSpectrumLogLikelihood[i][j]);
//				System.out.println(spectrumLogLikelihood[i][j] +"\t"+ scaledSpectrumLogLikelihood[i][j] +"\t"+ liS.getLogLikelihood());
//				}
			}	
			sumScaledSrpLogLikelihood[i] = liS.getSumScaledLikelihood();
			
			eachSrpLogLikelihood[i] = liS.getLogLikelihood();
//			System.out.println(i +"\t"+ eachSrpLogLikelihood[i]);
		}
		
//		double logLikelihood = liS.sumLogLikelihood(sumScaledSrpLogLikelihood);
		double totalLogLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return totalLogLikelihood;
	}



	
	private double calculateSrpLikelihoodSingle() {

//System.out.println("StartSingle");
		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int j = record.getSpectrumIndex(); 
		int k = record.getSingleIndex();//AllSiteIndexs()[0];
		double currentLogLikelihood = getStoredLogLikelihood();
		SpectraParameter spectra = spectrumModel.getSpectrum(j).getSpectra(k);
//		ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);

//		stateLikelihood.calculateStatesLogLikelihood2D(spectra, allStateLogLikelihood2D[0]);
//		stateLikelihood.calculateStoredStatesLogLikelihood2D(spectra, allStoredStateLogLikelihood2D[0]);

		stateLikelihood.calculateStatesLogLikelihood(spectra, 0, allStateLogLikelihood);
		stateLikelihood.calculateStoredStatesLogLikelihood(spectra, 0, allStoredStateLogLikelihood);

////		for (int i : mapToSrp) {
		int[] mapArray = mapToSrpArray[k];
		for (int i : mapArray){
//			String fullSrp = srpMap.getSrpFull(i);
//			int state = getStateAtK(fullSrp, k);
			int state = allSrpState2D[i][k];
			
//			currentLogLikelihood = updateLikelihoodAtIJK(i, j, state, allStateLogLikelihood2D[0],
//					allStoredStateLogLikelihood2D[0], currentLogLikelihood);
			if (state < STATE_COUNT) {
				currentLogLikelihood = updateLikelihoodAtIJK(i, j, state,
//						allStateLogLikelihood, allStoredStateLogLikelihood,
						currentLogLikelihood);
			}
		}

		return currentLogLikelihood;
	}
	
	
	private double calculateSrpLikelihoodMulti() {
		
		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();

		int[] siteIndexs = record.getAllSiteIndexs();
		int j= record.getSpectrumIndex(); 
		Spectrum spectrum = spectrumModel.getSpectrum(j);
		
//		stateLikelihood.calculateStatesLogLikelihood(spectra, allStateLogLikelihood2D[k]);
//		stateLikelihood.calculateStoredStatesLogLikelihood(spectra, allStoredStateLogLikelihood2D[k]);
		
		for (int k : siteIndexs) {
			SpectraParameter spectra = spectrum.getSpectra(k);
			int kOffset = k*STATE_COUNT;
			stateLikelihood.calculateStatesLogLikelihood(spectra, kOffset, allStateLogLikelihood);
			stateLikelihood.calculateStoredStatesLogLikelihood(spectra, kOffset, allStoredStateLogLikelihood);
		}
		
		int multihere;
		
		double currentLogLikelihood = getStoredLogLikelihood();
		if(multiType == MultiType.BitSet){
			recalculateBitSet(siteIndexs);
			
			int count = 0;
			for (int i = bitSet.nextSetBit(0); i >= 0; i = bitSet.nextSetBit(i+1)) {
				srpIndex[count++] = i;

				currentLogLikelihood = updateLikelihoodAtIJ_Local1DArray(i, j,
						siteIndexs, allStateLogLikelihood, allStoredStateLogLikelihood, currentLogLikelihood);

			}
			srpIndexCount = count;

			
		}
		else if(multiType==MultiType.Array){
			recalculateArray(siteIndexs);

			for (int i = 0; i < srpSwitch.length; i++) {
				if(srpSwitch[i]){

					currentLogLikelihood = updateLikelihoodAtIJ_Local1DArray(i, j, siteIndexs,
							allStateLogLikelihood, allStoredStateLogLikelihood, currentLogLikelihood);

				}
			}
	
		}
		else if(multiType==MultiType.Hash){
			recalculateHashSet(siteIndexs);

			for (int i : allSrpPos) {
				currentLogLikelihood = updateLikelihoodAtIJ_Local1DArray(i, j, siteIndexs,
						allStateLogLikelihood, allStoredStateLogLikelihood, currentLogLikelihood);
			}
		}
		else if(multiType==MultiType.All){
			for (int s = 0; s < siteIndexs.length; s++) {
				int k = siteIndexs[s];
				ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);

				for (int i : mapToSrp) {
//					String fullSrp = srpMap.getSrpFull(i);
//					int state = getStateAtK(fullSrp, k);
					int state = allSrpState2D[i][k];
					currentLogLikelihood = updateLikelihoodAtIJK(i, j, state,
							allStateLogLikelihood2D[s], allStoredStateLogLikelihood2D[s],
							currentLogLikelihood);
				}
			}
			
		}
		
		

		return currentLogLikelihood;

	}


	private void recalculateArray(int[] siteIndexs) {
		for (int s : siteIndexs) {
			for (int i : mapToSrpArray[s]){
				srpSwitch[i] = true;
			}
		}
	}

	private void recalculateHashSet(int[] siteIndexs) {
		allSrpPos.clear();
		for (int s : siteIndexs) {
			ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(s);
			allSrpPos.addAll(mapToSrp);
		}
	}

	private void recalculateBitSet(int[] siteIndexs) {
		bitSet.clear();
//		BitSet bitSet = new BitSet(srpCount);
		for (int s : siteIndexs) {
			BitSet tempSet = srpMap.getBitSet(s);
			bitSet.or(tempSet);
		}
	}

	private double calculateSrpLikelihoodColumn() {

		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int k = record.getSingleIndex();
		ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);
		int[] allSpectrumIndexs = record.getAllSpectrumIndexs();
		double currentLogLikelihood = getStoredLogLikelihood();
		
		for (int j : allSpectrumIndexs) {
			SpectraParameter spectra = spectrumModel.getSpectrum(j).getSpectra(k);
			
			stateLikelihood.calculateStatesLogLikelihood2D(spectra, allStateLogLikelihood2D[j]);
			stateLikelihood.calculateStoredStatesLogLikelihood2D(spectra, allStoredStateLogLikelihood2D[j]);
			int kOffset = j*STATE_COUNT;
			stateLikelihood.calculateStatesLogLikelihood(spectra, kOffset, allStateLogLikelihood);
			stateLikelihood.calculateStoredStatesLogLikelihood(spectra, kOffset, allStoredStateLogLikelihood);
		}
		for (int i : mapToSrp) {
			int state = allSrpState2D[i][k];
			for (int j : allSpectrumIndexs) {
//				currentLogLikelihood = updateLikelihoodAtIJK(i, j, state, allStateLogLikelihood2D[j],
//						allStoredStateLogLikelihood2D[j], currentLogLikelihood);
				if (state < STATE_COUNT) {
					currentLogLikelihood = updateLikelihoodAtIJK(i, j, j*STATE_COUNT+state,
							currentLogLikelihood);
				}
			}
		}
		return currentLogLikelihood;
	}
	
	private double calculateSrpLikelihoodSubColumn() {

		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int k = record.getSingleIndex();
		ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);
		int[] allSpectrumIndexs = record.getAllSpectrumIndexs();
		
		double currentLogLikelihood = getStoredLogLikelihood();

//		for (int j = 0; j < spectrumCount; j++) {
		for (int j : allSpectrumIndexs) {
			SpectraParameter spectra = spectrumModel.getSpectrum(j).getSpectra(k);
			
			stateLikelihood.calculateStatesLogLikelihood2D(spectra, allStateLogLikelihood2D[j]);
			stateLikelihood.calculateStoredStatesLogLikelihood2D(spectra, allStoredStateLogLikelihood2D[j]);
			int kOffset = j*spectrumCount;
			stateLikelihood.calculateStatesLogLikelihood(spectra, kOffset, allStateLogLikelihood);
			stateLikelihood.calculateStoredStatesLogLikelihood(spectra, kOffset, allStoredStateLogLikelihood);
	
		}
		
		for (int i : mapToSrp) {
//			String fullSrp = srpMap.getSrpFull(i);
//			int state = getStateAtK(fullSrp, k);
			int state = allSrpState2D[i][k];
//			for (int j = 0; j < spectrumCount; j++) {
			for (int j : allSpectrumIndexs) {
				currentLogLikelihood = updateLikelihoodAtIJK(i, j, state, allStateLogLikelihood2D[j],
						allStoredStateLogLikelihood2D[j], currentLogLikelihood);

			}


		}
		return currentLogLikelihood;
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
			ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				srpSwitch[i] = true;
			}
		}

		for (int k : siteIndexs) {
			
			SpectraParameter spectra0 = spectrum0.getSpectra(k);
			SpectraParameter spectra1 = spectrum1.getSpectra(k);
			int kOffset = k*STATE_COUNT;
			stateLikelihood.calculateStoredStatesLogLikelihood(spectra0, kOffset, allStateLogLikelihood);
			stateLikelihood.calculateStoredStatesLogLikelihood(spectra1, kOffset, allStoredStateLogLikelihood);
			
//			stateLikelihood.calculateStoredStatesLogLikelihood2D(spectra0, allStateLogLikelihood2D[s]);
//			stateLikelihood.calculateStoredStatesLogLikelihood2D(spectra1, allStoredStateLogLikelihood2D[s]);
			spectra0.setStateLikelihood(allStoredStateLogLikelihood, kOffset);
			spectra1.setStateLikelihood(allStateLogLikelihood, kOffset);
		
		}

		int multihere;
		double currentLogLikelihood = getStoredLogLikelihood();
		if(multiType == MultiType.BitSet){

			bitSet.clear();
//			BitSet bitSet = new BitSet(srpCount);
			for (int s : siteIndexs) {
				BitSet tempSet = srpMap.getBitSet(s);
				bitSet.or(tempSet);
			}
			
			int count = 0;
			for (int i = bitSet.nextSetBit(0); i >= 0; i = bitSet.nextSetBit(i+1)) {
				srpIndex[count] = i;
				count++;

				currentLogLikelihood = updateLikelihoodAtIJ_Local1DArray(i, j0,
						siteIndexs, allStoredStateLogLikelihood, allStateLogLikelihood, currentLogLikelihood);
				currentLogLikelihood = updateLikelihoodAtIJ_Local1DArray(i, j1,
						siteIndexs, allStateLogLikelihood, allStoredStateLogLikelihood, currentLogLikelihood);


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

					currentLogLikelihood = updateLikelihoodAtIJ_Local1DArray(i, j0,
							siteIndexs, allStoredStateLogLikelihood, allStateLogLikelihood, currentLogLikelihood);
					currentLogLikelihood = updateLikelihoodAtIJ_Local1DArray(i, j1,
							siteIndexs, allStateLogLikelihood, allStoredStateLogLikelihood, currentLogLikelihood);

				}
			}
	
		}
		

		
	
//		totalLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return currentLogLikelihood;

	}




	private double getStoredLogLikelihood() {
		return storedLogLikelihood;
	}
	
	private int getStateAtK(String fullSrp, int k) {
		char srpChar = fullSrp.charAt(k);//TODO: change to char[] at hight lv or at ShortReadMapping
		int state = DATA_TYPE.getState(srpChar);
		
		return state;
	}



	private double updateLikelihoodAtIJK(int i, int j, int state, 
//			double[] statesLogLikelihood, double[] storedStatesLogLikelihood, 
			double currentLogLikelihood) {
	
//		if(state<STATE_COUNT){
			double stateLn= allStateLogLikelihood[state];
			double storedStateLn = allStoredStateLogLikelihood[state];
			int offset = i*spectrumCount+j;
			if(storedStateLn != stateLn){
				currentLogLikelihood -= eachSrpLogLikelihood[i];
				sumScaledSrpLogLikelihood[i] -= scaledSpectrumLogLikelihood[offset];
		
				spectrumLogLikelihood[offset] -= storedStateLn; 
				spectrumLogLikelihood[offset] += stateLn;

				scaledSpectrumLogLikelihood[offset] = LikelihoodScaler.scale(spectrumLogLikelihood[offset], LOG_C);
				
				sumScaledSrpLogLikelihood[i] += scaledSpectrumLogLikelihood[offset];

				eachSrpLogLikelihood[i] = LikelihoodScaler.getLogLikelihood(sumScaledSrpLogLikelihood[i], LOG_C);
				currentLogLikelihood += eachSrpLogLikelihood[i];
			}

//			if(storedStateLn != stateLn){
//				currentLogLikelihood -= eachSrpLogLikelihood[i];
//				sumScaledSrpLogLikelihood[i] -= scaledSpectrumLogLikelihood2D[i][j];
//		
//				spectrumLogLikelihood2D[i][j] -= storedStateLn; 
//				spectrumLogLikelihood2D[i][j] += stateLn;
//		
//		//			scaledSpectrumLogLikelihood[i][j] = liS.scale(localSpectrum);
//				scaledSpectrumLogLikelihood2D[i][j] = LikelihoodScaler.scale(spectrumLogLikelihood2D[i][j], LOG_C);
//				
//				sumScaledSrpLogLikelihood[i] += scaledSpectrumLogLikelihood2D[i][j];
//		
//		//			eachSrpLogLikelihood[i] = liS.getLogLikelihood(sumScaledSrpLogLikelihood[i]);
//				eachSrpLogLikelihood[i] = LikelihoodScaler.getLogLikelihood(sumScaledSrpLogLikelihood[i], LOG_C);
//				currentLogLikelihood += eachSrpLogLikelihood[i];
//			}
//		}			
		return currentLogLikelihood;
	}

	private double updateLikelihoodAtIJK(int i, int j, int state, 
			double[] allStateLogLikelihood, double[] allStoredStateLogLikelihood, 
			double currentLogLikelihood) {
	
//		if(state<STATE_COUNT){
			double stateLn= allStateLogLikelihood[state];
			double storedStateLn = allStoredStateLogLikelihood[state];
			
			if(storedStateLn != stateLn){
				currentLogLikelihood -= eachSrpLogLikelihood[i];
				sumScaledSrpLogLikelihood[i] -= scaledSpectrumLogLikelihood2D[i][j];
		
				spectrumLogLikelihood2D[i][j] -= storedStateLn; 
				spectrumLogLikelihood2D[i][j] += stateLn;
		
		//			scaledSpectrumLogLikelihood[i][j] = liS.scale(localSpectrum);
				scaledSpectrumLogLikelihood2D[i][j] = LikelihoodScaler.scale(spectrumLogLikelihood2D[i][j], LOG_C);
				
				sumScaledSrpLogLikelihood[i] += scaledSpectrumLogLikelihood2D[i][j];
		
		//			eachSrpLogLikelihood[i] = liS.getLogLikelihood(sumScaledSrpLogLikelihood[i]);
				eachSrpLogLikelihood[i] = LikelihoodScaler.getLogLikelihood(sumScaledSrpLogLikelihood[i], LOG_C);
				currentLogLikelihood += eachSrpLogLikelihood[i];
		
			}
//		}			
		return currentLogLikelihood;
	}

	private double updateLikelihoodAtIJ_Local1DArray(int i, int j, int[] siteIndexs, 
			double[] allStateLogLikelihood, double[] allStoredStateLogLikelihood, double currentLogLikelihood) {

		int offsetIJ = i*spectrumCount+j;
		currentLogLikelihood -= eachSrpLogLikelihood[i];
		sumScaledSrpLogLikelihood[i] -= scaledSpectrumLogLikelihood[offsetIJ];
//		sumScaledSrpLogLikelihood[i] -= scaledSpectrumLogLikelihood2D[i][j];
		
//		double localSpectrum = spectrumLogLikelihood2D[i][j];
		double localSpectrum = spectrumLogLikelihood[offsetIJ];
		int[] thisSrpStata = allSrpState2D[i];
		for (int k : siteIndexs) {
			int state = thisSrpStata[k];
			if (state < STATE_COUNT) {
				int offset = k*STATE_COUNT+state;
//				double stateLn = allStateLogLikelihood[offset];
//				double storedStateLn = allStoredStateLogLikelihood[offset];
//				localSpectrum -= storedStateLn;
//				localSpectrum += stateLn;
//				
				localSpectrum -= allStoredStateLogLikelihood[offset];
				localSpectrum += allStateLogLikelihood[offset];
				

			}
		}
		
//		spectrumLogLikelihood2D[i][j] = localSpectrum;
//		scaledSpectrumLogLikelihood[i][j] = liS.scale(localSpectrum);
//		scaledSpectrumLogLikelihood2D[i][j] = LikelihoodScaler.scale(localSpectrum, LOG_C);
		
		spectrumLogLikelihood[offsetIJ] = localSpectrum;
		scaledSpectrumLogLikelihood[offsetIJ] = LikelihoodScaler.scale(localSpectrum, LOG_C);
		sumScaledSrpLogLikelihood[i] += scaledSpectrumLogLikelihood[offsetIJ];

//		eachSrpLogLikelihood[i] = liS.getLogLikelihood(sumScaledSrpLogLikelihood[i]);
		eachSrpLogLikelihood[i] = LikelihoodScaler.getLogLikelihood(sumScaledSrpLogLikelihood[i], LOG_C);
		currentLogLikelihood += eachSrpLogLikelihood[i];

		return currentLogLikelihood;
	}
	
	@Deprecated
	private double updateLikelihoodAtIJ(int i, int j, int[] siteIndexs, 
				double[][] allStateLogLikelihood, double[][] storedAllStateLogLikelihood, 
				double currentLogLikelihood) {

//		String fullSrp = srpMap.getSrpFull(i);
		
		currentLogLikelihood -= eachSrpLogLikelihood[i];
		liS.setsumScaledLikelihood(sumScaledSrpLogLikelihood[i]);
		liS.minus( scaledSpectrumLogLikelihood2D[i][j]);
		int[] thisSrpStata = allSrpState2D[i];
		for (int s = 0; s < siteIndexs.length; s++) {
			int k = siteIndexs[s];
//			int state = getStateAtK(fullSrp, k);
			int state = thisSrpStata[k];
			if(state<STATE_COUNT){
				double stateLn= allStateLogLikelihood[s][state];
				double storedStateLn = storedAllStateLogLikelihood[s][state];
				if(storedStateLn != stateLn){
					spectrumLogLikelihood2D[i][j] -= storedStateLn; 
					spectrumLogLikelihood2D[i][j] += stateLn;
				}
			}
		}

		scaledSpectrumLogLikelihood2D[i][j] = liS.scale(spectrumLogLikelihood2D[i][j]);
		liS.add(scaledSpectrumLogLikelihood2D[i][j]);
//		liS.addScaleLogProb(spectrumLogLikelihood[i][j]);
		sumScaledSrpLogLikelihood[i] = liS.getSumScaledLikelihood();
		eachSrpLogLikelihood[i] = liS.getLogLikelihood();
		currentLogLikelihood += eachSrpLogLikelihood[i];

		return currentLogLikelihood;
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

	@SuppressWarnings("rawtypes")
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
		ArrayList<Integer> mapToSrp;

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

		case DELTA_COLUMN:
		case SWAP_COLUMN:
		case SWAP_SUBCOLUMN:
			k= spectrumOperationRecord.getSingleIndex();
			mapToSrp = srpMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				storeI(i);
				for (j = 0; j < spectrumCount; j++) {
					storeIJ(i, j);
				}
			}

			break;
		case DELTA_SINGLE:
		case SWAP_SINGLE:
			j = spectrumOperationRecord.getSpectrumIndex();
			k = spectrumOperationRecord.getSingleIndex();//AllSiteIndexs()[0];
//			mapToSrp = srpMap.getMapToSrp(k);
//			spectrumModel.getSpectrum(j).getSpectra(k).storeState();
//			for (int i : mapToSrp) {
			for (int i : mapToSrpArray[k]){
				storeI(i);
				storeIJ(i, j);
//				storedSpectrumLogLikelihood[i][j] = spectrumLogLikelihood[i][j];
//				storedScaledSpectrumLogLikelihood[i][j] = scaledSpectrumLogLikelihood[i][j];
//				storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
//				storedSumSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i];

			}
			break;
			
		case DELTA_MULTI:
		case SWAP_MULTI:
			j = spectrumOperationRecord.getSpectrumIndex();
			int[] siteIndexs = spectrumOperationRecord.getAllSiteIndexs();
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
//				for (int s = 0; s < siteIndexs.length; s++) {
//					k = siteIndexs[s];
					mapToSrp = srpMap.getMapToSrp(kk);
					for (int i : mapToSrp) {
						storeI(i);
						storeIJ(i, j);
					}
				}
			}
						
			break;
		case RECOMBINATION:
			int[] twoSpectrums = spectrumOperationRecord.getRecombinationSpectrumIndex();
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
			throw new IllegalArgumentException("Unknown operation type: "+operation +"\tin"+ShortReadsSpectrumLikelihood.class.getSimpleName() );
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

	private void storeIJ(int i, int j) {
//		storedSpectrumLogLikelihood2D[i][j] = spectrumLogLikelihood2D[i][j];
//		storedScaledSpectrumLogLikelihood2D[i][j] = scaledSpectrumLogLikelihood2D[i][j];
		
		int offset = i*spectrumCount+j;
		storedSpectrumLogLikelihood[offset] = spectrumLogLikelihood[offset];
		storedScaledSpectrumLogLikelihood[offset] = scaledSpectrumLogLikelihood[offset];
	}

	private void storeI(int i) {
		storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
		storedSumSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i];
	}

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
		ArrayList<Integer> mapToSrp;
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
			
			restoreEverything();
			break;
		case DELTA_COLUMN:
		case SWAP_COLUMN:
		case SWAP_SUBCOLUMN:
			k = spectrumOperationRecord.getSingleIndex();
			mapToSrp = srpMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				restoreI(i);
				for (j = 0; j < spectrumCount; j++) {
					restoreIJ(i, j);
				}
			}

			break;
		case DELTA_SINGLE:
		case SWAP_SINGLE:
			j = spectrumOperationRecord.getSpectrumIndex();
			k = spectrumOperationRecord.getSingleIndex();//AllSiteIndexs()[0];
//			mapToSrp = srpMap.getMapToSrp(k);
//			spectrumModel.getSpectrum(j).getSpectra(k).restoreState();
//			for (int i : mapToSrp) {
			for (int i : mapToSrpArray[k]){
				restoreI(i);
				restoreIJ(i, j);
			}
			
			
		case DELTA_MULTI:
		case SWAP_MULTI:
			j = spectrumOperationRecord.getSpectrumIndex();
			int[] siteIndexs = spectrumOperationRecord.getAllSiteIndexs();
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
					mapToSrp = srpMap.getMapToSrp(kk);
					for (int i : mapToSrp) {
						restoreI(i);
						restoreIJ(i, j);
					}
				}

			}
			
			break;
		case RECOMBINATION:
			int[] twoSpectrumsIndex = spectrumOperationRecord.getRecombinationSpectrumIndex();

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

	private void restoreIJ(int i, int j) {
		int offset = i*spectrumCount+j;
		spectrumLogLikelihood[offset] = storedSpectrumLogLikelihood[offset];
		scaledSpectrumLogLikelihood[offset] = storedScaledSpectrumLogLikelihood[offset];

//		spectrumLogLikelihood2D[i][j] = storedSpectrumLogLikelihood2D[i][j];
//		scaledSpectrumLogLikelihood2D[i][j] = storedScaledSpectrumLogLikelihood2D[i][j];
	}

	private void restoreI(int i) {
		eachSrpLogLikelihood[i] = storedEachSrpLogLikelihood[i];
		sumScaledSrpLogLikelihood[i] = storedSumSrpLogLikelihood[i]; 
		
	}

	@Override
	protected void acceptState() {
		//Do nothing
	}
	private void storeEverything() {
	
			System.arraycopy(eachSrpLogLikelihood, 0, storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
			System.arraycopy(sumScaledSrpLogLikelihood, 0, storedSumSrpLogLikelihood, 0, sumScaledSrpLogLikelihood.length);
			System.arraycopy(spectrumLogLikelihood,0, storedSpectrumLogLikelihood, 0, spectrumLogLikelihood.length);
			System.arraycopy(scaledSpectrumLogLikelihood,0, storedScaledSpectrumLogLikelihood, 0, scaledSpectrumLogLikelihood.length);
			
			
			
			//		storedLogLikelihood = logLikelihood;
	//		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
	//		int spectrumIndex = record.getSpectrumIndex();
	//		int siteIndex = record.getSiteIndex();
	//		System.out.println(spectrumIndex +"\t"+ siteIndex);
			for (int i = 0; i < srpCount; i++) {
	//			System.err.println(i +"\t"+ allLogLikelihood[i][spectrumIndex][siteIndex] +"\t"+ storedAllLogLikelihood[i][spectrumIndex][siteIndex]);
	//			storedAllLogLikelihood[i][spectrumIndex][siteIndex] = allLogLikelihood[i][spectrumIndex][siteIndex];
				System.arraycopy(spectrumLogLikelihood2D[i],0, storedSpectrumLogLikelihood2D[i], 0, spectrumCount);
				System.arraycopy(scaledSpectrumLogLikelihood2D[i],0, storedScaledSpectrumLogLikelihood2D[i], 0, spectrumCount);
				for (int j = 0; j < spectrumCount; j++) {
	//				Spectrum spectrum = spectrumModel.getSpectrum(j);
	//				for (int s = 0; s < spectrumLength; s++) {
	//					spectrum.getSpectra(s).storeState();
	//				}
	//				for (int k = 0; k < spectrumLength; k++) {
	//					if(allLogLikelihood[i][j][k] != storedAllLogLikelihood[i][j][k]){
	//						System.out.println("DIFFLI:"+i +" "+j+" "+" "+k+
	//								" "+ allLogLikelihood[i][j][k] +" "+
	//								storedAllLogLikelihood[i][j][k]);
	//					}
	//				}
	//				System.arraycopy(allLogLikelihood[i][j], 0, storedAllLogLikelihood[i][j], 0, spectrumLength);
				}
				
			}
		}

	private void restoreEverything(){
		
		System.arraycopy(storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		System.arraycopy(storedSumSrpLogLikelihood, 0, sumScaledSrpLogLikelihood, 0, sumScaledSrpLogLikelihood.length);
		System.arraycopy(storedSpectrumLogLikelihood, 0, spectrumLogLikelihood, 0, storedSpectrumLogLikelihood.length);
		System.arraycopy(storedScaledSpectrumLogLikelihood, 0, scaledSpectrumLogLikelihood, 0, storedScaledSpectrumLogLikelihood.length);

		
		for (int i = 0; i < srpCount; i++) {
			
//			System.err.println(i +"\t"+ allLogLikelihood[i][spectrumIndex][siteIndex] +"\t"+ storedAllLogLikelihood[i][spectrumIndex][siteIndex]);
//			storedAllLogLikelihood[i][spectrumIndex][siteIndex] = allLogLikelihood[i][spectrumIndex][siteIndex];
			System.arraycopy(storedSpectrumLogLikelihood2D[i], 0, spectrumLogLikelihood2D[i], 0, spectrumCount);
			System.arraycopy(storedScaledSpectrumLogLikelihood2D[i], 0, scaledSpectrumLogLikelihood2D[i], 0, spectrumCount);
			for (int j = 0; j < spectrumCount; j++) {
//				for (int k = 0; k < spectrumLength; k++) {
//					if(allLogLikelihood[i][j][k] != storedAllLogLikelihood[i][j][k]){
//						System.out.println("DIFFLI:"+i +" "+j+" "+" "+k+
//								" "+ allLogLikelihood[i][j][k] +" "+
//								storedAllLogLikelihood[i][j][k]);
//					}
//				}
//				System.arraycopy(storedAllLogLikelihood[i][j], 0, allLogLikelihood[i][j], 0, spectrumLength);
			}
			
		}
	}
	public void setMultiType(MultiType type){
		this.multiType = type;
	}

	private void setDistType(DistType distType) {
//		this.distType = distType;
		try {

			switch (distType) {
			case flat:
				stateLikelihood = new ProbabilityStateLikelihood();
				break;
			case betaMean:
				stateLikelihood = new BetaMeanStateLikelihood();
				break;
			case betaMode:
				stateLikelihood = new BetaModeStateLikelihood();
				break;
			case gTest:
				stateLikelihood = new GTestStateLikelihood();
				break;
			case chisq:
				stateLikelihood = new ChisqStateLikelihood();
				break;
			default:
				throw new IllegalArgumentException("Invalid distType: " + distType);
			}
			likelihoodKnown = false;

		} catch (IllegalArgumentException e) {
			System.err.println("Invalid distribution type " + distType);
			e.printStackTrace();
			System.exit(-1);
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

	public double calculateSrpLikelihoodFullMaster() {


		double logLikelihood = 0;
		double spectrumLogLikelihood = 0;
		double stateLogLikelihood = 0;
		
		for (int i = 0; i < srpCount; i++) {

			String fullSrp = srpMap.getSrpFull(i);
			int start = srpMap.getSrpStart(i);
			int end = srpMap.getSrpEnd(i);
			
			liS.reset();
			for (int j = 0; j < spectrumCount; j++) {

				Spectrum spectrum = spectrumModel.getSpectrum(j);
				spectrumLogLikelihood = 0;
				for (int k = start; k < end; k++) {
					double[] frequencies = spectrum.getFrequenciesAt(k);
					int state = getStateAtK(fullSrp, k);
					if(state<STATE_COUNT){
//						stateLogLikelihood = caluclateStateLogLikelihood(frequencies[state]);
						stateLogLikelihood = stateLikelihood.caluclateStateLogLikelihood(frequencies[state]);
//						double likelihood = frequencies[state] * NOT_ERROR_RATE
//								+ (1 - frequencies[state]) * ERROR_RATE;
//						stateLogLikelihood = Math.log(likelihood);
					}
					else{
						stateLogLikelihood = MIN_LOG_LIKELIHOOD;
					}
					
					spectrumLogLikelihood += stateLogLikelihood;
					
				}
				liS.addLogProb(spectrumLogLikelihood);
				

			}	
			logLikelihood += liS.getLogLikelihood();
		}

//		double logLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		if(DEBUG){
			if(logLikelihood != this.logLikelihood){
				System.out.println(logLikelihood +"\t"+ this.logLikelihood +"\t"+ getStoredLogLikelihood());
			
			
			
				SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
//				int[] siteIndexs = record.getAllSiteIndexs();
				int j= record.getSpectrumIndex(); 
				Spectrum spectrum = spectrumModel.getSpectrum(j);
				
//				stateLikelihood.calculateStatesLogLikelihood(spectra, allStateLogLikelihood[k]);
//				stateLikelihood.calculateStoredStatesLogLikelihood(spectra, allStoredStateLogLikelihood[k]);
				
				for (int s = 0; s < spectrumLength; s++) {
					int k = s;
					SpectraParameter spectra = spectrum.getSpectra(k);
//					int kOffset = k*STATE_COUNT;
					double[] stateLn = spectra.getStateLikelihood();
//					stateLikelihood.calculateStatesLogLikelihood(spectra, 0, stateLn);
//					stateLikelihood.calculateStoredStatesLogLikelihood(spectra, kOffset, allStoredStateLogLikelihood);
//					System.out.println(Arrays.toString());
					
					double[] frequencies = spectrum.getFrequenciesAt(k);
					double[] stateLn2 = new double[4];
					for (int state = 0; state < 4; state++) {
						stateLn2[state] = stateLikelihood.caluclateStateLogLikelihood(frequencies[state]);
						if(stateLn2[state] != stateLn[state]){
							System.out.println(s);
							System.out.println(Arrays.toString(stateLn));
							System.out.println(Arrays.toString(stateLn2));
						}
					}
//					System.out.println(Arrays.toString(stateLn2));	
				}
//				System.exit(-1);
			
			
			
			
			
			
			
			}
		}
		return logLikelihood;
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


enum MultiType{
	Array,
	Hash,
	All, BitSet,;

};

public enum DistType{
	betaMean(0),
	betaMode(1),
	gTest(2), 
	chisq(3),
	flat(9),
	;
	int code;
	private DistType(int code) {
		this.code = code;
	}
//	DistType.valueOf(codeString);
}

//REMOVE after

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
				

 *//*

//REMOVE aftre
///////////////////////////////////////
@SuppressWarnings("unused")
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
				
				ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);
	//			System.out.println("Site: "+k +"\t"+ mapToSrp.size());
				for (int i : mapToSrp) {
					srpSwitch[i] = true;
				}
			}
			
			Set<Integer> allSrpPos = new HashSet<Integer>();
			allSrpPos.clear();
			for (int i = twoPositions[0]; i < twoPositions[1]; i++) {
				ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(i);
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


@Deprecated
private double calculateSrpLikelihoodSingle5() {


	SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
	int j = record.getSpectrumIndex(); 
	int k = record.getAllSiteIndexs()[0];

	Spectrum spectrum = spectrumModel.getSpectrum(j);
	ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);
	
	for (int i : mapToSrp) {
		String fullSrp = srpMap.getSrpFull(i);
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
			ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(s);
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
				String fullSrp = srpMap.getSrpFull(i);
			
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
	ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);

	
	
	for (int i : mapToSrp) {
		String fullSrp = srpMap.getSrpFull(i);
		
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
			ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				srpSwitch[i] = true;
			}
		}
		for (int i = 0; i < srpSwitch.length; i++) {
			if(srpSwitch[i]){
				
				String fullSrp = srpMap.getSrpFull(i);
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
			stateLikelihood.caluclateStateLogLikelihood(frequency);
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
*/
}
