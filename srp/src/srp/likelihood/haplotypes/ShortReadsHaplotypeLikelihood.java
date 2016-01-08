package srp.likelihood.haplotypes;

//import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.HashMap;

import javafx.scene.shape.CullFace;

import javax.swing.text.TabableView;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.ArithmeticUtils;

import com.sun.org.glassfish.external.statistics.Stats;









//import dr.math.BigDecimalUtils;
import dr.math.MathUtils;
import dr.math.Polynomial;
//import dr.math.Polynomial.BigDouble;
import dr.math.distributions.NormalDistribution;
//import srp.core.BigFunctions;
import srp.evolution.OperationType;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
//import srp.evolution.shortreads.AlignmentMapping;
import srp.evolution.shortreads.ShortReadMapping;
import srp.likelihood.AbstractShortReadsLikelihood;
import srp.likelihood.LikelihoodScaler;
import srp.likelihood.LikelihoodUtils;
//import srp.likelihood.spectrum.AbstractShortReadsSpectrumLikelihood;
import srp.likelihood.stateLikelihood.StateLikelihood;
//import java.util.BitSet;



public class ShortReadsHaplotypeLikelihood  extends AbstractShortReadsLikelihood {

	private static final long serialVersionUID = 7438385718398999755L;

	private static final boolean DEBUG = false;
	
//    public static final String SHORT_READ_LIKELIHOOD = "ShortReadHaplotypeLikelihood";
	public static final String NAME = "ShortReadHaplotypeLikelihood";

	public static final double decayRate = 0.49;
	private static final double LOG_DecayRate = Math.log(decayRate);
	private static final double LOG_ONE_MINUS_DecayRate = Math.log(1-decayRate);

	private static final double LOG_ONE_MINUS_RANDOM_RATE = Math.log(0.75);

	private static final double LOG_RANDOM_RATE = Math.log(0.25);
	

//	private final double MIN_LOG_LIKELIHOOD;

	

	
	protected HaplotypeModel alignmentModel;
	
	@Deprecated protected double[] allStateLogLikelihood;
	@Deprecated protected double[] allStoredStateLogLikelihood;

//	protected HaplotypeModel haplotypeModel;
	
	private HashMap<Integer, double[]> scaledLogBinomialDesnity;
//	private HashMap<Integer, BigDecimal[]> bigDecBinomialDesnity;
	private int[][] allDists;
	private int[][] storedAllDists;
	
	@Deprecated protected StateLikelihood stateLikelihood;

	public int globalCounter = 0;
	public int globalCounter2 = 0;
	public int globalCounter3 = 0;

	public int globalCounter4 = 0;

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
//		this.alignmentModel.storeModelState();
//		for (int j = 0; j < sequenceCount; j++) {
//			Haplotype haplotype = this.alignmentModel.getHaplotype(j);
//			haplotype.storeState();
//		}
		
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
		
		scaledLogBinomialDesnity = new HashMap<Integer, double[]>();
//		scaledLogBinomialDesnity = new double[]
		allDists = new int[srpCount][sequenceCount];
		storedAllDists = new int[srpCount][sequenceCount];

		double[] logFractorial = new double[500];//Change 500 later TODO:
		for (int i = 0; i < logFractorial.length; i++) {
			logFractorial[i] = ArithmeticUtils.factorialLog(i);
		}
		
		initScaledLikelihoodBinomial();
		
	}


	private void initScaledLikelihoodBinomial() {

		for (int s = 0; s < srpCount; s++) {
			int srLength = srpMap.getSrpLength(s);//ength();//srp.length();
			int srLengthPlusOne = srLength+1;
//			maxDist = Math.max(maxDist, srLength1);
			double[] binomD = new double[srLengthPlusOne];
			double[] logBinomD = new double[srLengthPlusOne];
			double[] scaledBinomD = new double[srLengthPlusOne];
//			double sum=0;
			
			if (!scaledLogBinomialDesnity.containsKey(srLength)) {
				for (int i = 0; i < logBinomD.length; i++) { // i is the number of mismatches (error) E^i 
				
					logBinomD[i] = ArithmeticUtils.binomialCoefficientLog(srLength, i)+i*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE;
					
//					double temp = ArithmeticUtils.binomialCoefficientLog(
//							srLength, i)
//							+ i
//							* LOG_RANDOM_RATE
//							+ (srLength - i)
//							* LOG_ONE_MINUS_RANDOM_RATE;
//					
//					logBinomD[i] = Math.log( Math.exp(logBinomD[i])*1/10 + Math.exp(temp)*9/10);
//					binomD[i] = Math.exp(logBinomD[i])/sequenceCount;
//					logBinomD[i] = Math.log(binomD[i]);
//					logBinomD[i] -= Math.log(10);
//					logBinomD[i] /= 7;
//					binomD[i] = Math.exp(logBinomD[i]);
//					logBinomD[i] = Math.log(binomD[i]);
//					System.out.println(i +"\t"+ logBinomD[i] +"\t"+ Math.log(binomD[i]) +"\t"+ (logBinomD[i] - Math.log(sequenceCount)) );
//				logBinomD[i] = Math.log(binomD[i]);
				
	
					scaledBinomD[i] = liS.scale(logBinomD[i]);
//					System.out.println(i +"\t"+ srLength +"\t"+ logBinomD.length +"\t"+ logBinomD[i] +"\t"+ Math.exp(logBinomD[i]) +"\t"+ scaledBinomD[i]);
//					sum += Math.exp(logBinomD[i]);
				}
			
			
//				System.out.println(s +"\t"+ srLength +"\t"+  logBinomD[0] +"\t"+ logBinomD[1] +"\t"+ logBinomD[2]);
//			System.exit(0);
				scaledLogBinomialDesnity.put(srLength, scaledBinomD);
			}
		}
		
	}


    
	protected double calculateSrpLikelihoodFull() {

//		System.out.println("calculateSrpLikelihood_Full");
		double logLikelihood = 0;
//		bigDecLikelihood = new BigDecimal("0");
		double overallSumDist = 0;
		for (int i = 0; i < srpCount; i++) {
			
			String srp = srpMap.getSrpFragment(i);//srpArray[i]
			int start = srpMap.getSrpStart(i);
			int end = srpMap.getSrpEnd(i);
			
			Integer srpLength = allSrpLengthInteger[i];
			double[] logPD = scaledLogBinomialDesnity.get(srpLength);

			liS.reset();
			sumScaledSrpLogLikelihood[i] = 0;
			for (int j = 0; j < sequenceCount; j++) {

				int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(j));
//				int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(j).getSequenceString());
				allDists[i][j]=dist;
				liS.add(logPD[dist]);
				overallSumDist += dist;
//				System.out.println(dist +"\t"+ logPD[dist]);
//				System.out.println(sumBigDecSrp[i]);
//				sumScaledSrpLogLikelihood[i] += logPD[allDists[i][j]];
//				storedAllLog2D[i][j] = logPD[dist];
//				liS.scaleLogProb(logPD[dist]);
			}	
//			System.out.println(Arrays.toString(allDists[i]));
//			System.out.println(liS.getSumScaledLikelihood());
//			System.out.println(liS.getLogLikelihood());
			sumScaledSrpLogLikelihood[i] = liS.getSumScaledLikelihood();
			eachSrpLogLikelihood[i] = liS.getLogLikelihood() ;
			
//			eachSrpLogLikelihood[i] = liS.getLogLikelihood();
//			bigDecEachSrpLogLikelihood[i] = BigFunctions.ln(sumBigDecSrp[i], 300);
//			eachSrpLogLikelihood[i] = bigDecEachSrpLogLikelihood[i].doubleValue();
//			System.out.println(eachSrpLogLikelihood[i] +"\t"+ sumBigDecSrp[i].toPlainString());
//			System.out.println("inFull: "+sumScaledSrpLogLikelihood[i] +"\t"+ eachSrpLogLikelihood[i]);
//			bigDecLikelihood = bigDecLikelihood.add(bigDecEachSrpLogLikelihood[i]);
			logLikelihood += eachSrpLogLikelihood[i];
		}
//		System.out.println("LogLikelihood: "+logLikelihood );
//		System.out.println(overallSumDist +"\t"+ srpCount +"\t");
		
		return logLikelihood;
	}

	private int calculateDeltaDist(char srpChar, char newChar, char oldChar) {
		int deltaDist = 0;

		if (srpChar == newChar) {
			deltaDist = -1;
		} else if (srpChar == oldChar) {
			deltaDist = 1;
		}
		// System.out.println("\tCalDelta: "+srpChar +"\t"+ newChar +"\t"+
		// oldChar +"\t"+ deltaDist);
		return deltaDist;
	}

	protected double calculateSrpLikelihoodSingle() {

		int j = operationRecord.getSpectrumIndex(); 
		int k = operationRecord.getSingleIndex();//AllSiteIndexs()[0];
		
		double currentLogLikelihood = getStoredLogLikelihood();
		double currentLogLikelihoodBackup = currentLogLikelihood;
		
		Haplotype haplotype = alignmentModel.getHaplotype(j);
		char oldChar = haplotype.getStoredChar(k);
		char newChar = haplotype.getChar(k);
		
//		System.out.println(currentLogLikelihood +"\t" +globalCounter);
globalCounter ++;
		if(newChar!= oldChar){ 
//		if(false){
			for (int i : mapToSrpArray[k]){
				char srpChar = allSrpChar2D[i][k];
				int deltaDist = calculateDeltaDist(srpChar, newChar, oldChar);
				
				if (deltaDist != 0) {
					
					allDists[i][j] += deltaDist;
					currentLogLikelihood = updateEachSrpAtIFull(i, currentLogLikelihood);
				}
				
					
			}
		}
		else{
			
		}
		if (currentLogLikelihoodBackup != currentLogLikelihood ) 
		{
			double delta = (currentLogLikelihoodBackup - currentLogLikelihood);
			double sum = 0;
			globalCounter2 ++;
			if(delta>0){
				globalCounter3++;
			}
			else{
				globalCounter4++;
			}
			for (int i : mapToSrpArray[k]){
				char srpChar = allSrpChar2D[i][k];
				int deltaDist = calculateDeltaDist(srpChar, newChar, oldChar);
				sum += deltaDist;
//				System.out.print(srpChar  +"\t");
			}
//			System.out.println("NOT Identical Likelihood " +"\t"+ oldChar +"\t"+ newChar +"\t"+ globalCounter2 +"/"+ globalCounter +"\t"+  delta );
//			System.out.println(sum +"\n"+ globalCounter2 +"\t"+ globalCounter3 +"\t"+ globalCounter4);

		}
//		currentLogLikelihood = bigDecLikelihood.doubleValue();
//		System.out.println(currentLogLikelihood);
		
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
			currentLogLikelihood = updateEachSrpAtIFull(i, currentLogLikelihood);
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
			throw new IllegalArgumentException("Deprecated");
		}
		

		
	
//		totalLikelihood = StatUtils.sum(eachSrpLogLikelihood);
		return currentLogLikelihood;

	}






	private double updateLikelihoodAtIJ(int i, int j, int[] siteIndexs, 
			Haplotype haplotype, double currentLogLikelihood) {

		
		updateAllDists(i, j, siteIndexs, haplotype);
//		int oldDist = allDists[i][j];
//		currentLogLikelihood = updateEachSrpAtI(i, currentLogLikelihood, oldDist, allDists[i][j]);
		currentLogLikelihood = updateEachSrpAtIFull(i, currentLogLikelihood);
		return currentLogLikelihood;
	}
	
	private double updateLikelihoodRecomb(int i, int j0, int j1, int[] siteIndexs, 
			Haplotype haplotype0, Haplotype haplotype1, double currentLogLikelihood) {

		updateAllDists(i, j0, siteIndexs, haplotype0);
//		currentLogLikelihood = updateEachSrpAtI(i, currentLogLikelihood);
		updateAllDists(i, j1, siteIndexs, haplotype1);
		currentLogLikelihood = updateEachSrpAtIFull(i, currentLogLikelihood);
		
//		updateLikelihoodAtIJ(i, j0, siteIndexs, haplotype0, currentLogLikelihood)
//		System.out.println(deltaCount +"\t"+ i);
		return currentLogLikelihood;
		
	}

	
	private double updateEachSrpAtIFull(int i, double currentLogLikelihood) {

		Integer srpLength = allSrpLengthInteger[i];
		double[] logPD = scaledLogBinomialDesnity.get(srpLength);
		int[] allD = allDists[i];
		
		double srpSum = 0;//~60 vs 90 for multi
		for (int h = 0; h < sequenceCount; h++) {
			srpSum += logPD[ allD[h] ];
		}
		
		currentLogLikelihood -= eachSrpLogLikelihood[i];
		eachSrpLogLikelihood[i] = LikelihoodScaler.getLogLikelihood(srpSum, LOG_C); //SLOW
		currentLogLikelihood += eachSrpLogLikelihood[i];
		

//		double xx = 0;
//		for (int hj = 0; hj < sequenceCount; hj++) {
////			xx *= Math.exp(logPD[js[hj]]);
//			xx+= Math.exp( LikelihoodScaler.getLogLikelihood(logPD[js[hj]], LOG_C)  );
//		}
//		
//		System.out.println(srpSum +"\t"+ eachSrpLogLikelihood[i] +"\t"+ Math.exp(eachSrpLogLikelihood[i]) +"\t"+ xx);
		
		return currentLogLikelihood;
	
	}
	
	/**
     * @deprecated
     * Lost precision in calculation.
     */
	@Deprecated
	private double updateEachSrpAtI(int i, double currentLogLikelihood, int oldDist, int newDist) {
		
		
		Integer srpLength = allSrpLengthInteger[i];
		double[] logPD = scaledLogBinomialDesnity.get(srpLength);
		
		double srpSum = 0;//~60 vs 90 for multi
		
		sumScaledSrpLogLikelihood[i] = storedSumScaledSrpLogLikelihood[i] - logPD[oldDist]+ logPD[newDist];
		srpSum = sumScaledSrpLogLikelihood[i];
		
		currentLogLikelihood -= eachSrpLogLikelihood[i];
		eachSrpLogLikelihood[i] = LikelihoodScaler.getLogLikelihood(srpSum, LOG_C);
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
			}
		}
//		int oldDist = storedAllDists[i][j];
//		int newDist = storedAllDists[i][j] + deltaDist;
		allDists[i][j] += deltaDist;
				
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
			System.out.println("StoreState in ShortReadsHaplotypeLikelihood:\t"+operation);
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
//			storeEverything();
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
			storeEverything();
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
				}throw new IllegalArgumentException("Deprecated");
			}
			else if(multiType==MultiType.Hash){
				for (int i : allSrpPos) {
					storeI(i);
					storeIJ(i, j);
				}throw new IllegalArgumentException("Deprecated");
			}
			else if(multiType==MultiType.All){
				for (int kk : siteIndexs) {
					for (int i : mapToSrpArray[kk]){
						storeI(i);
						storeIJ(i, j);
					}throw new IllegalArgumentException("Deprecated");
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
				}throw new IllegalArgumentException("Deprecated");
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

	
	private void storeEverything() {

		System.arraycopy(eachSrpLogLikelihood, 0, storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		System.arraycopy(sumScaledSrpLogLikelihood, 0, storedSumScaledSrpLogLikelihood, 0, sumScaledSrpLogLikelihood.length);
		
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
			System.out.println("RestoreState in ShortReadsHaplotypeLikelihood:\t"+operation);
		}
		switch (operation) {
		case NONE:
			
			break;
		case FULL:
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
//			restoreEverything();
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
			restoreEverything();
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
				}throw new IllegalArgumentException("Deprecated");
			}

			else if(multiType==MultiType.Hash){
				for (int i : allSrpPos) {
					restoreI(i);
					restoreIJ(i, j);
				}throw new IllegalArgumentException("Deprecated");
				
			}
			else if(multiType==MultiType.All){
				for (int kk : siteIndexs) {
					for (int i : mapToSrpArray[kk]){
						restoreI(i);
						restoreIJ(i, j);
					}
				}throw new IllegalArgumentException("Deprecated");

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
				}throw new IllegalArgumentException("Deprecated");
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
		
	}


	private void restoreIJ(int i, int j) {
		allDists[i][j] = storedAllDists[i][j];
		
	}


	private void restoreEverything() {
		
		System.arraycopy(storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		System.arraycopy(storedSumScaledSrpLogLikelihood, 0, sumScaledSrpLogLikelihood, 0, sumScaledSrpLogLikelihood.length);
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
			
			Integer srpLength = allSrpLengthInteger[i];
			double[] logPD = scaledLogBinomialDesnity.get(srpLength);

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
						
					for (int j = 0; j < sequenceCount; j++) {

						int dist = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(j).getSequenceString());
						int dist2 = LikelihoodUtils.Dist(start, end, srp, alignmentModel.getHaplotype(j));
						System.out.print(dist +"\t"+ dist2 +"\t"+ allDists[i][j]);
						if(dist!=dist2){
							System.out.print("dist!=dist2\t");
						}
						if(dist!=allDists[i][j]){
							System.out.print("dist!=allDists[i][j]\t");
						}
						if(dist2!=allDists[i][j]){
							System.out.print("dist2!=allDists[i][j]\t");
						}
						System.out.println();

					}
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



	private void updateAllDistsDebug(int i, int j, int[] siteIndexs, Haplotype haplotype) {
			
			char[] srpChars = allSrpChar2D[i];
			int deltaDist = 0;
			for (int k : siteIndexs) {
				char oldChar = haplotype.getStoredChar(k);
				char newChar = haplotype.getChar(k);
				if(newChar!= oldChar){ 
					deltaDist += calculateDeltaDist(srpChars[k], newChar, oldChar);
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


}
