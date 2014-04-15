package srp.likelihood.spectrum;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import javax.swing.text.TabableView;

import org.apache.commons.math3.stat.StatUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import srp.dr.evolution.datatype.ShortReads;
import srp.evolution.OperationType;
import srp.evolution.OperationRecord;
import srp.evolution.shortreads.ShortReadMapping;
import srp.evolution.spectrum.SpectraParameter;
import srp.evolution.spectrum.Spectrum;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import srp.likelihood.LikelihoodScaler;

import com.carrotsearch.hppc.BitSet;

import dr.evolution.datatype.DataType;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
//import java.util.BitSet;


public class ShortReadsSpectrumLikelihood  extends AbstractShortReadsSpectrumLikelihood {

	private static final long serialVersionUID = 7438385718398999755L;

	private static final boolean DEBUG = false;
	
//    public static final String SHORT_READ_LIKELIHOOD = "ShortReadSpectrumLikelihood";
	public static final String NAME = "ShortReadSpectrumLikelihood";

	private final double MIN_LOG_LIKELIHOOD;

//	private int srpCount;

	private LikelihoodScaler liS;

//	private int[][] allSrpState2D;
//
//	private int[] srpIndex;
//
//	private String[] srpArray;
//
//	private int srpIndexCount;

	protected double[] allStateLogLikelihood;
	protected double[] allStoredStateLogLikelihood;

	public ShortReadsSpectrumLikelihood(SpectrumAlignmentModel spectrumModel, ShortReadMapping srpMap){
		this(spectrumModel, srpMap, DistType.flat);
		
	}
	
	public ShortReadsSpectrumLikelihood(SpectrumAlignmentModel spectrumModel, ShortReadMapping srpMap, DistType distType){
		super(SHORT_READ_LIKELIHOOD, srpMap);
		this.spectrumModel = spectrumModel;

		operationRecord = spectrumModel.getOperationRecord();
//		multiType = MultiType.Array;
		multiType = MultiType.BitSet;

		setDistType(distType);
		
		MIN_LOG_LIKELIHOOD = stateLikelihood.caluclateStateLogLikelihood(SpectraParameter.MIN_FREQ);
		
		

		likelihoodKnown = false;
		
		addModel(this.spectrumModel);
		
		preprocessLikelihoodShortReadMap();
		calculateSrpLikelihoodFull();//TODO FIX this? shouldn't needed
		getLogLikelihood();
		
		
		storeEverything();
		
//		for (int i = 0; i < srpCount; i++) {
			for (int j = 0; j < sequenceCount; j++) {
				Spectrum spectrum = spectrumModel.getSpectrum(j);
				for (int s = 0; s < sequenceLength; s++) {
					spectrum.getSpectra(s).storeState();
				}
			}
//		}
	}
	

	private void preprocessLikelihoodShortReadMap() {
//		makeDirty();
		
		liS = new LikelihoodScaler(LOG_C);
		
		srpCount = srpMap.getSrpCount();
		sequenceCount = spectrumModel.getSpectrumCount();
		sequenceLength = spectrumModel.getSpectrumLength();
		
		this.srpSwitch = new boolean[srpCount];
		this.allSrpPos = new HashSet<Integer>();
		this.srpIndex = new int[srpCount];
		
		logLikelihood = Double.NEGATIVE_INFINITY;
		storedLogLikelihood = Double.NEGATIVE_INFINITY;

		spectrumLogLikelihood = new double[srpCount*sequenceCount];
		storedSpectrumLogLikelihood = new double[srpCount*sequenceCount];

		scaledSpectrumLogLikelihood = new double[srpCount*sequenceCount];
		storedScaledSpectrumLogLikelihood = new double[srpCount*sequenceCount];
		


		allStateLogLikelihood = new double[sequenceLength*STATE_COUNT];
		allStoredStateLogLikelihood = new double[sequenceLength*STATE_COUNT];
		
		mapToSrpArray = srpMap.getMapToSrpArray();
		
		String[] srpArray = srpMap.getSrpArray();
		
		allSrpState2D = new int[srpArray.length][sequenceLength];

		for (int i = 0; i < srpArray.length; i++) {
			String srp = srpArray[i];
			for (int j = 0; j < sequenceLength; j++) {
				allSrpState2D[i][j] = getStateAtK(srp, j);

			}
		}
		
		bitSet = new BitSet(srpCount);
	}


    
	protected double calculateSrpLikelihoodFull() {

//		System.out.println("calculateSrpLikelihoodFull");

		double[][] allStateLogLikelihoodFull2D = new double[sequenceCount][sequenceLength*STATE_COUNT];
		
		for (int j = 0; j < sequenceCount; j++) {
			Spectrum spectrum = spectrumModel.getSpectrum(j);
			for (int k = 0; k < sequenceLength; k++) {
				SpectraParameter spectra = spectrum.getSpectra(k);
				int kOffset = k*STATE_COUNT;
				stateLikelihood.calculateStatesLogLikelihood(spectra, kOffset, allStateLogLikelihood);
			}
			System.arraycopy(allStateLogLikelihood, 0  , allStateLogLikelihoodFull2D[j], 0, allStateLogLikelihood.length );
		}

		double stateLogLikelihood;
		for (int i = 0; i < srpCount; i++) {

//			String fullSrp = srpMap.getSrpFull(i);
			int start = srpMap.getSrpStart(i);
			int end = srpMap.getSrpEnd(i);
			
			liS.reset();
			for (int j = 0; j < sequenceCount; j++) {
				int offset = i*sequenceCount+j;
//				Spectrum spectrum = spectrumModel.getSpectrum(j);

				spectrumLogLikelihood[offset] = 0;
				for (int k = start; k < end; k++) {
//					int state = getStateAtK(fullSrp, k);
					int state = allSrpState2D[i][k];
					
					if(state<STATE_COUNT){
						int kOffset = k * STATE_COUNT+state;
						 stateLogLikelihood = allStateLogLikelihoodFull2D[j][kOffset];
					}
					else{
						stateLogLikelihood = MIN_LOG_LIKELIHOOD;
					}
//					System.out.println(stateLogLikelihood);
//					allLogLikelihood[i][j][k] = stateLogLikelihood;
					spectrumLogLikelihood[offset] += stateLogLikelihood;
					
				}
				
				scaledSpectrumLogLikelihood[offset] = liS.scale(spectrumLogLikelihood[offset]);
				
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



	
	protected double calculateSrpLikelihoodSingle() {

//System.out.println("StartSingle");
//		OperationRecord record = spectrumModel.getOperationRecord();
		int j = operationRecord.getSpectrumIndex(); 
		int k = operationRecord.getSingleIndex();//AllSiteIndexs()[0];
		double currentLogLikelihood = getStoredLogLikelihood();
		SpectraParameter spectra = spectrumModel.getSpectrum(j).getSpectra(k);
//		ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);

//		stateLikelihood.calculateStatesLogLikelihood2D(spectra, allStateLogLikelihood2D[0]);
//		stateLikelihood.calculateStoredStatesLogLikelihood2D(spectra, allStoredStateLogLikelihood2D[0]);

		stateLikelihood.calculateStatesLogLikelihood(spectra, 0, allStateLogLikelihood);
		stateLikelihood.calculateStoredStatesLogLikelihood(spectra, 0, allStoredStateLogLikelihood);

		for (int i : mapToSrpArray[k]){
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
	
	
	protected double calculateSrpLikelihoodMulti() {
		
//		OperationRecord record = spectrumModel.getOperationRecord();

		int[] siteIndexs = operationRecord.getAllSiteIndexs();
		int j= operationRecord.getSpectrumIndex(); 
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
		
		
		

		return currentLogLikelihood;

	}


	protected double calculateSrpLikelihoodColumn() {

//		OperationRecord record = spectrumModel.getOperationRecord();
		int k = operationRecord.getSingleIndex();
		ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);
		int[] allSpectrumIndexs = operationRecord.getAllSpectrumIndexs();
		double currentLogLikelihood = getStoredLogLikelihood();
		
		for (int j : allSpectrumIndexs) {
			SpectraParameter spectra = spectrumModel.getSpectrum(j).getSpectra(k);
			
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
	
	protected double calculateSrpLikelihoodSubColumn() {

		OperationRecord record = spectrumModel.getOperationRecord();
		int k = record.getSingleIndex();
		ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(k);
		int[] allSpectrumIndexs = record.getAllSpectrumIndexs();
		
		double currentLogLikelihood = getStoredLogLikelihood();

//		for (int j = 0; j < spectrumCount; j++) {
		for (int j : allSpectrumIndexs) {
			SpectraParameter spectra = spectrumModel.getSpectrum(j).getSpectra(k);
			
			int kOffset = j*sequenceCount;
			stateLikelihood.calculateStatesLogLikelihood(spectra, kOffset, allStateLogLikelihood);
			stateLikelihood.calculateStoredStatesLogLikelihood(spectra, kOffset, allStoredStateLogLikelihood);
	
		}
		
		for (int i : mapToSrp) {
//			String fullSrp = srpMap.getSrpFull(i);
//			int state = getStateAtK(fullSrp, k);
			int state = allSrpState2D[i][k];
//			for (int j = 0; j < spectrumCount; j++) {
			for (int j : allSpectrumIndexs) {
				if (state < STATE_COUNT) {
					currentLogLikelihood = updateLikelihoodAtIJK(i, j, j*STATE_COUNT+state,
							currentLogLikelihood);
				}
			}
			

		}
		return currentLogLikelihood;
	}

	

	protected double calculateSrpLikelihoodRecombination() {
		
		OperationRecord record = spectrumModel.getOperationRecord();
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




	



	private double updateLikelihoodAtIJK(int i, int j, int state, 
//			double[] statesLogLikelihood, double[] storedStatesLogLikelihood, 
			double currentLogLikelihood) {
	
//		if(state<STATE_COUNT){
			double stateLn= allStateLogLikelihood[state];
			double storedStateLn = allStoredStateLogLikelihood[state];
			int offset = i*sequenceCount+j;
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

	private double updateLikelihoodAtIJ_Local1DArray(int i, int j, int[] siteIndexs, 
			double[] allStateLogLikelihood, double[] allStoredStateLogLikelihood, double currentLogLikelihood) {

		int offsetIJ = i*sequenceCount+j;
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
	

	@Override
	protected void storeState() {
//long time1 = System.currentTimeMillis();

//		System.arraycopy(eachSrpLogLikelihood, 0, storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		storedLogLikelihood = logLikelihood;
//		storeEverything();
		OperationRecord spectrumOperationRecord = spectrumModel.getOperationRecord();
		OperationType operation = spectrumOperationRecord.getOperation();
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

		case COLUMN:
		
		case SWAP_SUBCOLUMN:
			k= spectrumOperationRecord.getSingleIndex();
			mapToSrp = srpMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				storeI(i);
				for (j = 0; j < sequenceCount; j++) {
					storeIJ(i, j);
				}
			}

			break;
		case SINGLE:
		
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
			
		case MULTI:
		
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


	@Override
	protected void restoreState() {
//		long time1 = System.currentTimeMillis();
		
//		System.err.println("SR likelihood restore: " + logLikelihood +"\t"+ storedLogLikelihood);
		logLikelihood = storedLogLikelihood;
//		restoreEverything();
//		System.arraycopy(storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		
		OperationRecord spectrumOperationRecord = spectrumModel.getOperationRecord();
		OperationType operation = spectrumOperationRecord.getOperation();
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
		case COLUMN:
		
		case SWAP_SUBCOLUMN:
			k = spectrumOperationRecord.getSingleIndex();
			mapToSrp = srpMap.getMapToSrp(k);
			for (int i : mapToSrp) {
				restoreI(i);
				for (j = 0; j < sequenceCount; j++) {
					restoreIJ(i, j);
				}
			}

			break;
		
		case SINGLE:
			j = spectrumOperationRecord.getSpectrumIndex();
			k = spectrumOperationRecord.getSingleIndex();//AllSiteIndexs()[0];
//			mapToSrp = srpMap.getMapToSrp(k);
//			spectrumModel.getSpectrum(j).getSpectra(k).restoreState();
//			for (int i : mapToSrp) {
			for (int i : mapToSrpArray[k]){
				restoreI(i);
				restoreIJ(i, j);
			}
			
			
		
		case MULTI:
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


	public double[] unittestMethodGetEachLikelihood() {
		double[] copyOfValues = new double[eachSrpLogLikelihood.length];
		System.arraycopy(eachSrpLogLikelihood, 0, copyOfValues, 0,
				copyOfValues.length);
		return copyOfValues;
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
			for (int j = 0; j < sequenceCount; j++) {

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
			
			
			
//				OperationRecord record = spectrumModel.getOperationRecord();
//				int[] siteIndexs = record.getAllSiteIndexs();
				int j= operationRecord.getSpectrumIndex(); 
				Spectrum spectrum = spectrumModel.getSpectrum(j);
				
//				stateLikelihood.calculateStatesLogLikelihood(spectra, allStateLogLikelihood[k]);
//				stateLikelihood.calculateStoredStatesLogLikelihood(spectra, allStoredStateLogLikelihood[k]);
				
				for (int s = 0; s < sequenceLength; s++) {
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


}
