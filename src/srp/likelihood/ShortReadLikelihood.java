package srp.likelihood;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.ArithmeticUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.Operation;
import srp.haplotypes.ShortRead;
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
//	public static final int[] NULL_SWAPINFO = HaplotypeModel.NULL_SWAPINFO;
	
//	public static final int POISSON = 0;
//	public static final int BINOMIAL = 0;
	
//	private double[] allFactorialLog;
//	private int maxDist;
	

	private int hapLength;
	private int haplotypeCount;
	private int srpCount;
	
    
	
	private double logLikelihood;
    private double storedLogLikelihood;
    protected boolean likelihoodKnown = false;
    
    private double[] eachSrpLikelihood;
    private double[] storedEachSrpLikelihood;
	
    private int[][] allDists;
	private int[][] storedAllDists;
	
	

	private HashMap<Integer, double[]> logBinomialDesnity;
	private HashMap<Integer, double[]> scaledLogBinomialDesnity;
	
	
	private ArrayList<char[]> shortReadChars;
	
	@Deprecated
	private int[] counter;
//	private HashMap<Integer, double[]> logPoissonDesnity;
//	SwapInfo swapInfo = new SwapInfo();
	private AlignmentMapping aMap;
	private HaplotypeModel haplotypeModel;
		

	private LikelihoodScaler liS;
	private double EVALUATION_TEST_THRESHOLD = 1e-8;
	
	
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


	public ShortReadLikelihood(HaplotypeModel haplotypeModel){
		this(SHORT_READ_LIKELIHOOD);
		this.haplotypeModel = haplotypeModel;
		this.aMap = this.haplotypeModel.getAlignmentMapping();
		
		
		preprocessLikelihoodAlignmentMap();
		calculateSrpLikelihoodFullUseCounter();
		
		addModel(this.haplotypeModel);
	}
	

	private void preprocessLikelihoodAlignmentMap() {
		makeDirty();
		
		liS = new LikelihoodScaler(LOG_C);
		
		srpCount = aMap.getSrpCount();
		haplotypeCount = haplotypeModel.getHaplotypeCount();
		allDists = new int[srpCount][haplotypeCount];
		storedAllDists = new int[srpCount][haplotypeCount];

		eachSrpLikelihood = new double[srpCount];
		storedEachSrpLikelihood = new double[srpCount];

		
		shortReadChars = new ArrayList<char[]>();

//		allFactorialLog = new double[hapLength+1];
//		for (int i = 0; i < allFactorialLog.length; i++) {
//			allFactorialLog[i] = ArithmeticUtils.factorialLog(i);
//		}
		
		

		logBinomialDesnity = new HashMap<Integer, double[]>();
		scaledLogBinomialDesnity = new HashMap<Integer, double[]>();
		
		int maxDist=0;
		for (int s = 0; s < srpCount; s++) {
			String srp = aMap.getSrpFragment(s);
			shortReadChars.add(srp.toCharArray());
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
		Operation op = haplotypeModel.getOperation();
		return calculateLogLikelihoodSelect(op);
	}

	public double calculateLogLikelihoodSelect(Operation op){

		double logLikelihood = Double.NEGATIVE_INFINITY;

		switch (op) {
			case NONE:
				logLikelihood = calculateSrpLikelihoodFull2();
				break;
			case SWAPMULTI:
				logLikelihood = calculateSrpLikelihoodMultiBasesSwap2();
//				logLikelihood = calculateSrpLikelihoodFull();
//				logLikelihood = storedLogLikelihood;
				break;
			case SWAPSINGLE:
				logLikelihood = calculateSrpLikelihoodSingleBaseSwap2();
//				logLikelihood = calculateSrpLikelihoodFull2();
//				logLikelihood = storedLogLikelihood;
				break;
//			case UNIFORMSWAPBASE:
//				logLikelihood = calculateSrpLikelihoodSingleBaseSwap();
//				break;
			case RECOMB:
				logLikelihood = calculateSrpLikelihoodFull2();
//				logLikelihood = storedLogLikelihood;
				break;
			case SWAPSECTION:
				logLikelihood = calculateSrpLikelihoodFull2();
//				logLikelihood = storedLogLikelihood;
				break;
	
			default:
				throw new IllegalArgumentException("Unknown operation type: "+op);
	
			}
//	    double logLikelihood = calculateShoreReadLikelihood4();
//	    double logLikelihood = calculateShoreReadLikelihoodBinomialModel2();
	    
//	    timeTrial();
//		storeState();
		return logLikelihood;
	}


	@Override
	public Model getModel() {
		return this;
		
	}



	@Override
	public void makeDirty() {
        likelihoodKnown = false;
		
	}

	
	private double calculateSrpLikelihoodFull2() {

		liS.reset();
		for (int i = 0; i < srpCount; i++) {

			char[] srpChar = shortReadChars.get(i);
			int start = aMap.getSrpStart(i);
			int end = aMap.getSrpEnd(i);
			double[] logPD = scaledLogBinomialDesnity.get(aMap.getSrpLength(i));
//			double[] logPD = logBinomialDesnity.get(aMap.getSrpLength(i));
			
			for (int j = 0; j < haplotypeCount; j++) {

				int dist = LikelihoodUtils.Dist(start, end, srpChar, haplotypeModel.getAlignedSequenceString(j));
				allDists[i][j]=dist;
				liS.addScaledLogProb(logPD[dist]);
//				liS.scaleLogProb(logPD[dist]);
			}	
			
			eachSrpLikelihood[i] = liS.getLogLikelihood();
			liS.reset();

		}	

		double logLikelihood = StatUtils.sum(eachSrpLikelihood);
//		double logLikelihood2 = calculateSrpLikelihoodFull();
//		if(logLikelihood != logLikelihood2){
//			System.out.println(logLikelihood +"\t"+ logLikelihood2 +"\t"+ (logLikelihood-logLikelihood2));
//		}
//		
		
		return logLikelihood;
	}


	private double calculateSrpLikelihoodSingleBaseSwap2() {
		
		int[] swapInfo = haplotypeModel.getSwapInfo().getSwapInfoSWAPBASE();
		int hapIndex = swapInfo[0];
		int swapPos = swapInfo[1];
		int newChar = swapInfo[2];
		int oldChar = swapInfo[3];
		
		String newHaplotypeStirng = haplotypeModel.getAlignedSequenceString(hapIndex);
		int hapChar = haplotypeModel.getHaplotypeCharAt(hapIndex, swapPos);
		ArrayList<Integer> mapPos = aMap.getMapToSrp(swapPos);
		
		boolean isHapEqualNew = hapChar==newChar;
		
		for (int srpIndex : mapPos) {
//			calLikeliSrpHap(srpIndex, hapIndex, newHaplotypeStirng);
			calLikeliSrpHapSingle(srpIndex, hapIndex, swapPos, newChar, oldChar, isHapEqualNew);
			
/*
			int start = aMap.getSrpStart(srpIndex);
			int end = aMap.getSrpEnd(srpIndex);
			String srpString = aMap.getSrpFragment(srpIndex);
//			double[] logPD = logBinomialDesnity.get(aMap.getSrpLength(i));
			double[] logPD = scaledLogBinomialDesnity.get(aMap.getSrpLength(srpIndex));
			
			int dist = LikelihoodUtils.Dist(start, end, srpString, newHaplotypeStirng );
			allDists[srpIndex][hapIndex] = dist;
			
			for (int j = 0; j < allDists[srpIndex].length; j++) {
//				liS.scaleLogProb(logPD[allDists[i][j]]);
				liS.addScaledLogProb(logPD[allDists[srpIndex][j]]);
			}
			eachSrpLikelihood[srpIndex] = liS.getLogLikelihood();
			liS.reset();
*/		
		}

		double logLikelihood = StatUtils.sum(eachSrpLikelihood);
		double logLikelihood2 = calculateSrpLikelihoodFullUseCounter();

		if(Math.abs( logLikelihood - logLikelihood2) > EVALUATION_TEST_THRESHOLD  ){
			System.out.println("SS: "+logLikelihood +"\t"+ logLikelihood2 +"\t"+ (logLikelihood-logLikelihood2));
		}
		return logLikelihood;
	}


	private double calculateSrpLikelihoodMultiBasesSwap2() {
		
		liS.reset();
		Deque<int[]> swapMulti = haplotypeModel.getSwapInfo().getSwapInfoSWAPMULTI();
		int[] swapInfo = swapMulti.peekLast();
		
//		HashSet<Integer> mapPos = new HashSet<Integer>();
//		for (Iterator<int[]> iterator = swapMulti.descendingIterator(); iterator
//				.hasNext();) {
//
//			nextSwap = iterator.next();
//
//			int swapPos = nextSwap[1];
//			ArrayList<Integer> mapToSrp = aMap.getMapToSrp(swapPos);
//			mapPos.addAll(mapToSrp);
//		}
		
//		int[] swapInfo = haplotypeModel.getSwapInfo().getSwapInfoSWAPBASE();
		int hapIndex = swapInfo[0];
		int swapPos = swapInfo[1];
		int newChar = swapInfo[2];
		int oldChar = swapInfo[3];
		
		String newHaplotypeStirng = haplotypeModel.getAlignedSequenceString(hapIndex);
		int hapChar = haplotypeModel.getHaplotypeCharAt(hapIndex, swapPos);
		
		boolean isHapEqualNew = hapChar==newChar;

		
//		ArrayList<Integer> mapPos = aMap.getMapToSrp(swapPos);
//		for (int i : mapPos) {
		for (int srpIndex = 0; srpIndex < srpCount; srpIndex++) {
				

				calLikeliSrpHap(srpIndex, hapIndex, newHaplotypeStirng);
//				calLikeliSrpHapSingle(srpIndex, hapIndex, swapPos, newChar, oldChar, isHapEqualNew);
				
			
		}

		double logLikelihood = StatUtils.sum(eachSrpLikelihood);
		double logLikelihood2 = calculateSrpLikelihoodFullUseCounter();
		if(Math.abs( logLikelihood - logLikelihood2) > EVALUATION_TEST_THRESHOLD  ){
			System.out.println("SS: "+logLikelihood +"\t"+ logLikelihood2 +"\t"+ (logLikelihood-logLikelihood2));
		}
		return logLikelihood;
	}


	private void calLikeliSrpHapSingle(int srpIndex, int hapIndex, int swapPos, 
			int newChar, int oldChar, boolean isHapEqualNew){

//		int start = srp.getStart();
//		int end = srp.getEnd();
//		String srpString = srp.getFragmentSrp();
//		int dist = LikelihoodUtils.Dist(start, end, srpString, newHaplotypeStirng );

		ShortRead srp = aMap.getShortRead(srpIndex);
		int srpChar = srp.getFullSrpCharAt(swapPos);
		int deltaDist = 0;
		if(newChar!= oldChar && isHapEqualNew){
			if (srpChar==newChar){
				deltaDist = -1;
			}
			else if(srpChar==oldChar){
				deltaDist = 1;
			}
		}
		if (deltaDist!= 0){
			double[] logPD = scaledLogBinomialDesnity.get(srp.getLength());

			int newDist = storedAllDists[srpIndex][hapIndex] + deltaDist;
	
			allDists[srpIndex][hapIndex] = newDist;
	
			liS.reset();		
			for (int j = 0; j < haplotypeCount ; j++) {
				liS.addScaledLogProb(logPD[allDists[srpIndex][j]]);
			}
			eachSrpLikelihood[srpIndex] = liS.getLogLikelihood();
		}
		
		
	}
	
	
	private void calLikeliSrpHap(int srpIndex, int hapIndex, String newHaplotypeStirng){
		
		
		ShortRead srp= aMap.getShortRead(srpIndex);
		int start = srp.getStart();
		int end = srp.getEnd();
		String srpString = srp.getFragmentSrp();
		double[] logPD = scaledLogBinomialDesnity.get(srp.getLength());

		int dist = LikelihoodUtils.Dist(start, end, srpString, newHaplotypeStirng );
		allDists[srpIndex][hapIndex] = dist;

		liS.reset();		
		for (int j = 0; j < haplotypeCount ; j++) {
			liS.addScaledLogProb(logPD[allDists[srpIndex][j]]);
		}
		eachSrpLikelihood[srpIndex] = liS.getLogLikelihood();
		

		
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
		// TODO Auto-generated method stub
		
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
		
//		System.arraycopy(storedEachLikelihood, 0, eachLikelihood, 0, eachLikelihood.length);
		
	}

	@Deprecated
		private double calculateSrpLikelihoodFullUseCounter() {
	
			double logLikelihood = 0;
			
			liS.reset();
			
			for (int i = 0; i < srpCount; i++) {
	
				Arrays.fill(counter, 0);
				char[] srpChar = shortReadChars.get(i);
				int start = aMap.getSrpStart(i);
				int end = aMap.getSrpEnd(i);
	//	System.out.println("=========\ni_start_end\t"+i +"\t"+ start +"\t"+ end);
	//			for (char[] hapCharArray : haplotypesChars) {
				for (int j = 0; j < haplotypeCount; j++) {
					
	//				char[] hapCharArray = haplotypesChars[j]; 
	
					int dists = LikelihoodUtils.Dist(start, end, srpChar, haplotypeModel.getAlignedSequenceString(j));
					counter[dists]++;
	
				}	
				
				double[] logPD = logBinomialDesnity.get(aMap.getSrpLength(i));
				for (int j = 0; j < counter.length; j++) {
					if(counter[j]!=0){
	//	System.out.println(j +" outOf "+(end-start) +"\t"+ logPD[j] +"\t"+ counter[j]);
						liS.scaleLogProb(logPD[j], counter[j]);
					}
				}
				
	//	System.out.println("===============");
				eachSrpLikelihood[i] = liS.getLogLikelihood();
	//			tempEachLikelihood[i] =eachLikelihood[i] ;
				logLikelihood += eachSrpLikelihood[i];
				liS.reset();
			}	
	//		storeState();
			
			
			return logLikelihood;
		}

	@Deprecated
		private double calculateSrpLikelihoodSingleBaseSwapOld() {
	//		System.out.println("Partial");
	//		System.out
	//		.println("\n==========\nEQUAL Likelihood:" +"\t"+ Arrays.equals(eachLikelihood, tempEachLikelihood));
	//		for (int i = 0; i < tempEachLikelihood.length; i++) {
	//			System.out.println(i +"\t"+ tempEachLikelihood[i] +"\t"+ eachLikelihood[i]);
	//		}
			
	
			double logLikelihood = 0;
			liS.reset();
			
	
			int swapPos = haplotypeModel.getSwapInfo().getSwapInfoSWAPBASE()[1];
			ArrayList<Integer> mapPos = aMap.getMapToSrp(swapPos);
			
			
			for (int i : mapPos) {
	
				Arrays.fill(counter, 0);
				char[] srpChar = shortReadChars.get(i);
				int start = aMap.getSrpStart(i);
				int end = aMap.getSrpEnd(i);
				String srpString = aMap.getSrpFragment(i);
				
				for (int j = 0; j < haplotypeCount; j++) {
	//				int dist = LikelihoodUtils.Dist(start, end, srpChar, haplotypeModel.getAlignedSequenceString(j));
	//				int dist = LikelihoodUtils.Dist(start, end, srpChar, haplotypeModel.getAlignedSequenceString(j).toCharArray());
					int dist = LikelihoodUtils.Dist(start, end, srpString, haplotypeModel.getAlignedSequenceString(j));
	//				int dist = j;
					counter[dist]++;
				}	
				
				double[] logPD = logBinomialDesnity.get(aMap.getSrpLength(i));
				for (int dist = 0; dist < logPD.length; dist++) {
					if(counter[dist]!=0){
						liS.scaleLogProb(logPD[dist], counter[dist]);
					}
					
				}
				eachSrpLikelihood[i] = liS.getLogLikelihood();
				liS.reset();
			}
	
			logLikelihood = StatUtils.sum(eachSrpLikelihood);
	
	//		double tempLL = StatUtils.sum(tempL);
	//		
	//		if (tempLL != logLikelihood)	{
	//			System.out.println("EXIT");
	//			
	//			System.exit(-1);
	//			System.out.println((tempLL == logLikelihood) + "\t" + tempLL + "\t"
	//					+ logLikelihood + "\t" + Arrays.toString(swapInfo) + "\t"+ mapPos.toString());
	//			for (int i = 0; i < tempL.length; i++) {
	//				if (tempL[i]!= tempEachLikelihood[i]){
	//					System.out.println(i +"\t"+ tempL[i] +"\t"+ tempEachLikelihood[i]);
	//				}
	//				
	//			}
	//			
	//
	////			for (int i = 0; i < haplotypesChars.length; i++) {
	////				System.out.println(Arrays.toString(haplotypesChars[i]));
	////			}
	//			System.out.println(mapPos.toString());
	//////			System.out.println(Arrays.toString(swapInfo)+"\t"+ mapPos.size());
	////			for (int i = 0; i < haplotypesChars.length; i++) {
	////				System.out.println(Arrays.toString(haplotypesChars[i]));
	////			}
	//			
	//		}
			
	//		storeState();
			return logLikelihood;
		}

	@Deprecated
	private double calculateSrpLikelihoodMultiBasesSwapOld() {
			
	
			
			liS.reset();
			Deque<int[]> swapMulti = haplotypeModel.getSwapInfo().getSwapInfoSWAPMULTI();
	
	//		ArrayList<Integer> mapPos = new ArrayList<Integer>();
			HashSet<Integer> mapPos = new HashSet<Integer>();
			for (Iterator<int[]> iterator = swapMulti.descendingIterator(); iterator
					.hasNext();) {
	
				int swapPos = iterator.next()[1];
				ArrayList<Integer> mapToSrp = aMap.getMapToSrp(swapPos);
				
				mapPos.addAll(mapToSrp);
	//			System.out.println(mapToSrp.toString());
			}
	//		System.out.println(mapPos.toString());
	//		System.exit(-1);
	//		int swapPos = 0;
			
	//		ArrayList<Integer> mapPos = aMap.getMapToSrp(swapPos);
			
			
			for (int i : mapPos) {
						
				Arrays.fill(counter, 0);
				char[] srpChar = shortReadChars.get(i);
				int start = aMap.getSrpStart(i);
				int end = aMap.getSrpEnd(i);
	
				for (int j = 0; j < haplotypeCount; j++) {
	//				int dists = LikelihoodUtils.Dist(start, end, srpChar, haplotypesChars[j]);
					int dists = LikelihoodUtils.Dist(start, end, srpChar, haplotypeModel.getAlignedSequenceString(j));
					counter[dists]++;
				}	
				
				double[] logPD = logBinomialDesnity.get(aMap.getSrpLength(i));
				for (int j = 0; j < counter.length; j++) {
					if(counter[j]!=0){
						liS.scaleLogProb(logPD[j], counter[j]);
					}
				}
	
				eachSrpLikelihood[i] = liS.getLogLikelihood();
				liS.reset();
			}
	
			double logLikelihood = StatUtils.sum(eachSrpLikelihood);
	//		double tempLL = StatUtils.sum(tempL);
	//		
	//		if (tempLL != logLikelihood)	{
	//			System.out.println("EXIT");
	//			
	//			System.exit(-1);
	//			System.out.println((tempLL == logLikelihood) + "\t" + tempLL + "\t"
	//					+ logLikelihood + "\t" + Arrays.toString(swapInfo) + "\t"+ mapPos.toString());
	//			for (int i = 0; i < tempL.length; i++) {
	//				if (tempL[i]!= tempEachLikelihood[i]){
	//					System.out.println(i +"\t"+ tempL[i] +"\t"+ tempEachLikelihood[i]);
	//				}
	//				
	//			}
	//			
	//
	////			for (int i = 0; i < haplotypesChars.length; i++) {
	////				System.out.println(Arrays.toString(haplotypesChars[i]));
	////			}
	//			System.out.println(mapPos.toString());
	//////			System.out.println(Arrays.toString(swapInfo)+"\t"+ mapPos.size());
	////			for (int i = 0; i < haplotypesChars.length; i++) {
	////				System.out.println(Arrays.toString(haplotypesChars[i]));
	////			}
	//			
	//		}
			
	//		storeState();
			return logLikelihood;
		
			
		}

}