package srp.likelihood.haplotypes;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.ArithmeticUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import srp.haplotypes.old.OldHaplotype;
import srp.haplotypes.old.OldHaplotypeModel;
import srp.haplotypes.old.OldHapOperation;
import srp.haplotypes.old.OldHapSwapInfo;
import srp.likelihood.LikelihoodScaler;
import srp.shortreads.AlignmentMapping;
import srp.shortreads.ShortRead;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;


@SuppressWarnings({ "unused", "serial" })
public class ShortReadLikelihood extends AbstractModelLikelihood {

	
    public static final String SHORT_READ_LIKELIHOOD = "ShortReadLikelihood";
	public static final String NAME = SHORT_READ_LIKELIHOOD;
	public static final double ERROR_RATE = 0.0107;
//	public static final double ONE_MINUS_ERROR_RATE = 1-ERROR_RATE;
	public static final double LOG_ERROR_RATE = Math.log(ERROR_RATE);
	public static final double LOG_ONE_MINUS_ERROR_RATE = Math.log(1-ERROR_RATE);
	public static final double C = 1e-200;
	public static final double LOG_C = Math.log(C);
	
	private static final double EVALUATION_TEST_THRESHOLD = 1e-8;
	//	public static final int[] NULL_SWAPINFO = HaplotypeModel.NULL_SWAPINFO;
	protected boolean likelihoodKnown = false;
	
	private int haplotypeLength;
	private int haplotypeCount;
	private int srpCount;

	private double logLikelihood;
	private double storedLogLikelihood;
	

	private double[] eachSrpLikelihood;
	private double[] storedEachSrpLikelihood;

	private int[][] allDists;
	private int[][] storedAllDists;

	private int[] allDelta;


	private HashMap<Integer, double[]> logBinomialDesnity;
	private HashMap<Integer, double[]> scaledLogBinomialDesnity;

	private AlignmentMapping aMap;
	private OldHaplotypeModel haplotypeModel;
	private LikelihoodScaler liS;
	private OldHapSwapInfo swapInfo;
	
	@Deprecated
	private int[] counter;
	private OldHapOperation operation;
	
	
	@Override
	public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

    public ShortReadLikelihood(String name) {

        super(SHORT_READ_LIKELIHOOD);
        likelihoodKnown = false;
        setId(SHORT_READ_LIKELIHOOD);
        
//
//        this.trialsParameter = trialsParameter;
//        this.proportionParameter = proportionParameter;
//        addVariable(trialsParameter);
//        addVariable(proportionParameter);
//        this.counts = counts;

    }


	public ShortReadLikelihood(OldHaplotypeModel haplotypeModel){
		this(SHORT_READ_LIKELIHOOD);
		this.haplotypeModel = haplotypeModel;
		this.aMap = this.haplotypeModel.getAlignmentMapping();
		this.swapInfo = this.haplotypeModel.getSwapInfo();
		operation = OldHapOperation.NONE;
		preprocessLikelihoodAlignmentMap();
		calculateSrpLikelihoodFull();
//		calculateSrpLikelihoodFullUseCounter();
		
		addModel(this.haplotypeModel);
	}
	

	private void preprocessLikelihoodAlignmentMap() {
		makeDirty();
		
		liS = new LikelihoodScaler(LOG_C);
		
		srpCount = aMap.getSrpCount();
		haplotypeCount = haplotypeModel.getHaplotypeCount();
		haplotypeLength = haplotypeModel.getHaplotypeLength();
		
		logLikelihood = Double.NEGATIVE_INFINITY;
		storedLogLikelihood = Double.NEGATIVE_INFINITY;
		
		allDists = new int[srpCount][haplotypeCount];
		storedAllDists = new int[srpCount][haplotypeCount];

		eachSrpLikelihood = new double[srpCount];
		storedEachSrpLikelihood = new double[srpCount];

		logBinomialDesnity = new HashMap<Integer, double[]>();
		scaledLogBinomialDesnity = new HashMap<Integer, double[]>();
		
		allDelta = new int[srpCount];

		
		int maxDist=0;
		for (int s = 0; s < srpCount; s++) {
			String srp = aMap.getSrpFragment(s);

			int srLength = srp.length();
//			char[] srCharArray = reads.toCharArray();
//			double plambda = srLength * ERROR_RATE;
//			double logPlambda = Math.log(plambda);

			int srLength1 = srLength+1;
			maxDist = Math.max(maxDist, srLength1);
			double[] logBinomD = new double[srLength1];
			double[] scaledBinomD = new double[srLength1];
			for (int i = 0; i < logBinomD.length; i++) {

				logBinomD[i] = ArithmeticUtils.binomialCoefficientLog(srLength, i)+i*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE;
//				logBinomD[i] = i*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE;
				scaledBinomD[i] = liS.scale(logBinomD[i]); 
			}
//			System.out.println(Arrays.toString(logBinomD));
			logBinomialDesnity.put(srLength, logBinomD);
			scaledLogBinomialDesnity.put(srLength, scaledBinomD);
		}
		
		
		counter = new int[maxDist];
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
		OldHapSwapInfo swapInfo = haplotypeModel.getSwapInfo();
		operation = swapInfo.getOperation();

		double logLikelihood = Double.NEGATIVE_INFINITY;

		switch (operation) {
			case NONE:
				logLikelihood = calculateSrpLikelihoodFull();
				break;
			case SWAPMULTI:
				logLikelihood = calculateSrpLikelihoodMultiBasesSwap();
				break;

			case SWAPSINGLE:
				logLikelihood = calculateSrpLikelihoodSingleBaseSwap(swapInfo);
				break;
//			case UNIFORMSWAPBASE:
//				logLikelihood = calculateSrpLikelihoodSingleBaseSwap();
//				break;
			case SWAPCOLUMN:
				logLikelihood = calculateSrpLikelihoodSwapColumn();
				break;
			case SWAPSECTION:
				logLikelihood = calculateSrpLikelihoodSwapSection();
				break;
			case PASS:
				logLikelihood = storedLogLikelihood;
				break;
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
		for (int i = 0; i < srpCount; i++) {

			String srp = aMap.getSrpFragment(i);
			int start = aMap.getSrpStart(i);
			int end = aMap.getSrpEnd(i);
			
			double[] logPD = scaledLogBinomialDesnity.get(aMap.getSrpLength(i));

			liS.reset();
			for (int j = 0; j < haplotypeCount; j++) {

				int dist = LikelihoodUtils.Dist(start, end, srp, haplotypeModel.getAlignedSequenceString(j));
				allDists[i][j]=dist;
				liS.add(logPD[dist]);
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
		
		return logLikelihood;
	}


	private double calculateSrpLikelihoodSingleBaseSwap(OldHapSwapInfo swapInfo) {
		
		int[] swapBaseRecord = swapInfo.getSwapInfoSWAPBASE();
		int hapIndex = swapBaseRecord[0];
		int swapPos = swapBaseRecord[1];
		int newChar = swapBaseRecord[2];
		int oldChar = swapBaseRecord[3];
		
		ArrayList<Integer> mapPos = aMap.getMapToSrp(swapPos);
		
		for (int srpIndex : mapPos) {

			calLikeliSrpHapSingle(srpIndex, hapIndex, swapPos, newChar, oldChar);
		}	
//		}

		double logLikelihood = StatUtils.sum(eachSrpLikelihood);

//		double logLikelihood2 = calculateSrpLikelihoodFullUseCounter();
//		if(Math.abs( logLikelihood - logLikelihood2) > EVALUATION_TEST_THRESHOLD  ){
//			System.out.println("LL_Single: "+logLikelihood +"\t"+ logLikelihood2 +"\t"+ (logLikelihood-logLikelihood2));
//		}
		
		return logLikelihood;
	}


	private double calculateSrpLikelihoodMultiBasesSwap() {
		

		int hapIndex = haplotypeModel.getSwapInfo().getHapIndex();
		int[][] swapMulti = haplotypeModel.getSwapInfo().getSwapInfoSWAPMULTI();
		int[] allNewChars = swapMulti[0];
		int[] allOldChars = swapMulti[1];
		
		Arrays.fill(allDelta, 0);

		for (int swapPos = 0; swapPos < haplotypeLength ; swapPos++) {
			int oldChar = allOldChars[swapPos];
			if(oldChar>0){
				int newChar = allNewChars[swapPos];
				ArrayList<Integer> mapPos = aMap.getMapToSrp(swapPos);
				for (int srpIndex : mapPos) {
	
					allDelta[srpIndex] += calculateDeltaDist(aMap.getShortReadCharAt(srpIndex, swapPos), newChar, oldChar);//, isHapEqualNew);
				}
			}
		}
		
		for (int srpIndex = 0; srpIndex < srpCount; srpIndex++) {
			if ( allDelta[srpIndex] != 0){
				
				int newDist = storedAllDists[srpIndex][hapIndex] + allDelta[srpIndex];
				allDists[srpIndex][hapIndex] = newDist;
		
				liS.reset();
				double[] logPD = scaledLogBinomialDesnity.get(aMap.getSrpLength(srpIndex));
				for (int j = 0; j < haplotypeCount ; j++) {
					liS.add(logPD[allDists[srpIndex][j]]);
				}
				eachSrpLikelihood[srpIndex] = liS.getLogLikelihood();
			}
		}
		

		double logLikelihood = StatUtils.sum(eachSrpLikelihood);

//		double logLikelihood2 = calculateSrpLikelihoodFullUseCounter();
//		if(Math.abs( logLikelihood - logLikelihood2) > EVALUATION_TEST_THRESHOLD  ){
//			System.out.println("LL_Multi: "+logLikelihood +"\t"+ logLikelihood2 +"\t"+ (logLikelihood-logLikelihood2));
//		}
		
		return logLikelihood;
	}
	


	private double calculateSrpLikelihoodSwapColumn() {
		
		int[][] swapColumn = haplotypeModel.getSwapInfo().getSwapInfoSWAPCOLUMN();
		int[] posChar = swapColumn[0];
		int[] allOldChars = swapColumn[1];
		
		int swapPos = posChar[0];
		int newChar = posChar[1];
		ArrayList<Integer> mapPos = aMap.getMapToSrp(swapPos);

		for (int srpIndex : mapPos) {

			ShortRead srp = aMap.getShortRead(srpIndex);
			int srpChar = srp.getFullSrpCharAt(swapPos);
			double[] logPD = scaledLogBinomialDesnity.get(srp.getLengthInteger());

//			boolean srpEqNew = srpChar==newChar;
			//TODO: is this faster? on average case?

			liS.reset();
			for (int hapIndex = 0; hapIndex < haplotypeCount; hapIndex++) {
				
				int oldChar = allOldChars[hapIndex];
				int deltaDist = calculateDeltaDist(srpChar, newChar, oldChar);
//				int deltaDist = calculateDeltaDist(srpEqNew, srpChar, newChar, oldChar);
			
				if (deltaDist!= 0){
					int newDist = storedAllDists[srpIndex][hapIndex] + deltaDist;
					allDists[srpIndex][hapIndex] = newDist;

				}
				liS.add(logPD[allDists[srpIndex][hapIndex]]);
			}
			eachSrpLikelihood[srpIndex] = liS.getLogLikelihood();

		}	
		
		double logLikelihood = StatUtils.sum(eachSrpLikelihood);

//		double logLikelihood2 = calculateSrpLikelihoodFull();
//		if(Math.abs( logLikelihood - logLikelihood2) > EVALUATION_TEST_THRESHOLD  ){
//			System.out.println("LL_Column: "+logLikelihood +"\t"+ logLikelihood2 +"\t"+ (logLikelihood-logLikelihood2));
//		}
		

		return logLikelihood;
	}



	private double calculateSrpLikelihoodSwapSection() {
		

		int[] swapHapRecord = haplotypeModel.getSwapInfo().getSwapInfoSWAPSECTION();
		int hapIndex1 = swapHapRecord[0];
		int hapIndex2 = swapHapRecord[1];
		
		Arrays.fill(allDelta, 0);

		OldHaplotype haplotype1 = haplotypeModel.getHaplotype(hapIndex1);
		OldHaplotype haplotype2 = haplotypeModel.getHaplotype(hapIndex2);
		for (int swapPos = swapHapRecord[2]; swapPos < swapHapRecord[3]; swapPos++) {
			int char1 = haplotype1.getChar(swapPos);
			int char2 = haplotype2.getChar(swapPos);
			
			ArrayList<Integer> mapPos = aMap.getMapToSrp(swapPos);
			for (int srpIndex : mapPos) {
				char srpChar = aMap.getShortReadCharAt(srpIndex, swapPos);
				allDelta[srpIndex] += calculateDeltaDist(srpChar, char1, char2);
				// delta of hapIndex1 == reverse of hapIndex2
			}
		}
		
		for (int srpIndex = 0; srpIndex < srpCount; srpIndex++) {
			if ( allDelta[srpIndex] != 0 ){
				
				int newDist = storedAllDists[srpIndex][hapIndex1] + allDelta[srpIndex];
				allDists[srpIndex][hapIndex1] = newDist;
				
				int newDist2 = storedAllDists[srpIndex][hapIndex2] - allDelta[srpIndex];
				allDists[srpIndex][hapIndex2] = newDist2;

				liS.reset();
				double[] logPD = scaledLogBinomialDesnity.get(aMap.getSrpLength(srpIndex));
				for (int j = 0; j < haplotypeCount ; j++) {
					liS.add(logPD[allDists[srpIndex][j]]);
				}
				eachSrpLikelihood[srpIndex] = liS.getLogLikelihood();
			}
			
		}
		

		double logLikelihood = StatUtils.sum(eachSrpLikelihood);

//		double logLikelihood2 = calculateSrpLikelihoodFullUseCounter();
//		if(Math.abs( logLikelihood - logLikelihood2) > EVALUATION_TEST_THRESHOLD  ){
//			System.out.println("LL_Section: "+logLikelihood +"\t"+ logLikelihood2 +"\t"+ (logLikelihood-logLikelihood2));
//		}
//		
		return logLikelihood;
	}
	
	private void calLikeliSrpHapSingle(int srpIndex, int hapIndex, int swapPos, 
			int newChar, int oldChar){

//		int start = srp.getStart();
//		int end = srp.getEnd();
//		String srpString = srp.getFragmentSrp();
//		int dist = LikelihoodUtils.Dist(start, end, srpString, newHaplotypeStirng );

		ShortRead srp = aMap.getShortRead(srpIndex);
		int srpChar = srp.getFullSrpCharAt(swapPos);
//		System.out.println(srpIndex +"\t"+  hapIndex+"\t"+  swapPos+"\t"+  srpChar +"\t"+newChar+"\t"+  oldChar );
		int deltaDist = calculateDeltaDist(srpChar, newChar, oldChar);

		if (deltaDist!= 0){
			double[] logPD = scaledLogBinomialDesnity.get(srp.getLength());

			int newDist = storedAllDists[srpIndex][hapIndex] + deltaDist;
//			System.out.println(storedAllDists[srpIndex][hapIndex] +"\t"+ deltaDist +"\t"+ newChar +"\t"+ oldChar +"\t"+ srpChar);
			allDists[srpIndex][hapIndex] = newDist;
	
			liS.reset();		
			for (int j = 0; j < haplotypeCount ; j++) {
				liS.add(logPD[allDists[srpIndex][j]]);
			}
			eachSrpLikelihood[srpIndex] = liS.getLogLikelihood();
		}

		
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
	private static int calculateDeltaDist(boolean srpEqNew, int srpChar, int newChar, int oldChar){//, boolean isHapEqualNew){
		
		int deltaDist = 0;
	
		if(newChar!= oldChar){ // if(newChar!= oldChar && isHapEqualNew)
			if (srpEqNew){
				deltaDist = -1;
			}
			else if(srpChar==oldChar){
				deltaDist = 1;
			}
		}
		return deltaDist;
		
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == haplotypeModel) {
            // treeModel has changed so recalculate the intervals
//            eventsKnown = false;
        }

        likelihoodKnown = false;
		
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {

		
	}

	@Deprecated
	public void restoreStatePublicTest(){
		restoreState();
	}

	@Override
	protected void storeState() {
//		System.err.println("SR likelihood store: " + logLikelihood);
		System.arraycopy(eachSrpLikelihood, 0, storedEachSrpLikelihood, 0, eachSrpLikelihood.length);
		storedLogLikelihood = logLikelihood;
//		System.arraycopy(storedEachLikelihood, 0, eachLikelihood, 0, eachLikelihood.length);

		for (int i = 0; i < allDists.length; i++) {
		    System.arraycopy(allDists[i], 0, storedAllDists[i], 0, allDists[0].length);
		}
		
		
	}
	@Override
	protected void restoreState() {
//		System.err.println("SR likelihood restore: " + storedLogLikelihood);
		logLikelihood = storedLogLikelihood;
		System.arraycopy(storedEachSrpLikelihood, 0, eachSrpLikelihood, 0, eachSrpLikelihood.length);
		

		for (int i = 0; i < storedAllDists.length; i++) {
		    System.arraycopy(storedAllDists[i], 0, allDists[i], 0, storedAllDists[0].length);
		}
		
		
	}

	@Override
	protected void acceptState() {
		//Do nothing
//		System.out.println("Accept");
	}

	public OldHapOperation getOperation(){
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


	@Deprecated
	private double calculateSrpLikelihoodFullUseCounter() {

		for (int i = 0; i < srpCount; i++) {

			Arrays.fill(counter, 0);

			int start = aMap.getSrpStart(i);
			int end = aMap.getSrpEnd(i);
			String srpString = aMap.getSrpFragment(i);

			for (int j = 0; j < haplotypeCount; j++) {

				int dists = LikelihoodUtils.Dist(start, end, srpString,
						haplotypeModel.getAlignedSequenceString(j));
				counter[dists]++;

			}
			liS.reset();
			double[] logPD = logBinomialDesnity.get(aMap.getSrpLength(i));
			for (int j = 0; j < counter.length; j++) {
				if (counter[j] != 0) {
					liS.addLogProbMulti(logPD[j], counter[j]);
				}
			}

			eachSrpLikelihood[i] = liS.getLogLikelihood();

		}
		logLikelihood = StatUtils.sum(eachSrpLikelihood);
		return logLikelihood;
	}

}