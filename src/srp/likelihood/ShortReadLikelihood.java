package srp.likelihood;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import javax.sound.midi.SysexMessage;


import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.ArithmeticUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentModel;
import srp.haplotypes.HaplotypeModel;

import com.google.common.cache.AbstractCache.StatsCounter;


import dr.evolution.alignment.Alignment;
import dr.evolution.sequence.Sequences;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;


@SuppressWarnings({ "unused", "serial" })
public class ShortReadLikelihood extends AbstractModelLikelihood {

	
    public static final String SHORT_READ_LIKELIHOOD = "ShortReadLikelihood";
	public static final String NAME = SHORT_READ_LIKELIHOOD;
	public static final double ERROR_RATE = 0.0107;
	public static final double ONE_MINUS_ERROR_RATE = 1-ERROR_RATE;
	public static final double LOG_ERROR_RATE = Math.log(ERROR_RATE);
	public static final double LOG_ONE_MINUS_ERROR_RATE = Math.log(1-ERROR_RATE);
	public static final double C = 1e-200;
	public static final double LOG_C = Math.log(C);
	public static final int[] NULL_SWAPINFO = HaplotypeModel.NULL_SWAPINFO;
	
//	public static final int POISSON = 0;
	public static final int BINOMIAL = 0;
	
	private double logLikelihood;
    private double storedLogLikelihood;
    protected boolean likelihoodKnown = false;
	
    private double[] eachLikelihood;
    private double[] storedEachLikelihood;
	
	

	private int maxDist;
	private double[] allFactorialLog;
	
	private char[][] oldHaplotypesChars;
	private char[][] haplotypesChars;
	private ArrayList<char[]> shortReadChars;
	
	
//	private HashMap<Integer, double[]> logPoissonDesnity;
	private HashMap<Integer, double[]> logBinomialDesnity;
	private Parameter tempValue;
	private AlignmentModel alignmentModel;
	private int hapLength;

	private int[] swapInfo;
	
	private AlignmentMapping aMap;
	
	
	@Override
	public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

    public ShortReadLikelihood(String name) {

        super(SHORT_READ_LIKELIHOOD);
        likelihoodKnown = false;
        setId(SHORT_READ_LIKELIHOOD);
        swapInfo = new int[4];
//
//        this.trialsParameter = trialsParameter;
//        this.proportionParameter = proportionParameter;
//        addVariable(trialsParameter);
//        addVariable(proportionParameter);
//        this.counts = counts;

    }

	public ShortReadLikelihood(AlignmentMapping aMap, HaplotypeModel haplotypeModel){
		this(SHORT_READ_LIKELIHOOD);
		
		

//		ArrayList<String> haplotypes = new ArrayList<>();
//		for (int i = 0; i < alignment.getSequenceCount(); i++) {
//			haplotypes.add(alignment.getAlignedSequenceString(i));
//		}
//
//
//		setHaplotypes(haplotypes);

		makeDirty();
		this.aMap = aMap;
		updateHaplotypes(haplotypeModel);
//		this.haplotypesChars = alignment.getCharMatrix();
		preprocessLikelihoodAlignmentMap();
		calculateShoreReadLikelihoodBinomalFull();
		
		addModel(haplotypeModel);
	}
	
	public void updateHaplotypes(HaplotypeModel aMatrix){
//		updatehaplotypesChars(aMatrix, null);
//	}
//	public void updatehaplotypesChars(Haplotypes aMatrix, int[] swapInfo){
		
		haplotypesChars = aMatrix.getCharMatrix();
		swapInfo = aMatrix.getSwapInfo();
//		for (int i = 0; i < haplotypesChars.length; i++) {
//			boolean isSame = Arrays.equals(oldHaplotypesChars[i], haplotypesChars[i]);
//			if (!isSame){
//				for (int j = 0; j < array.length; j++) {
//					
//				}
//			}
//					
//		}
		
		makeDirty();
	}
	

	private void preprocessLikelihoodAlignmentMap() {
		makeDirty();
		shortReadChars = new ArrayList<char[]>();

		allFactorialLog = new double[hapLength+1];
		for (int i = 0; i < allFactorialLog.length; i++) {
			allFactorialLog[i] = ArithmeticUtils.factorialLog(i);
		}
		
		maxDist=0;

		logBinomialDesnity = new HashMap<Integer, double[]>();
//		for (String reads : shortRead) {
		eachLikelihood = new double[aMap.getSrpCount()];
		storedEachLikelihood = new double[aMap.getSrpCount()];
		
		
		for (int s = 0; s < aMap.getSrpCount(); s++) {
			String srp = aMap.getSrpFragment(s);
			shortReadChars.add(srp.toCharArray());
			int srLength = srp.length();
//			char[] srCharArray = reads.toCharArray();
//			double plambda = srLength * ERROR_RATE;
//			double logPlambda = Math.log(plambda);

			int srLength1 = srLength+1;
			maxDist = Math.max(maxDist, srLength1);
			double[] logBinomD = new double[srLength1];
			for (int i = 0; i < logBinomD.length; i++) {

				logBinomD[i] = ArithmeticUtils.binomialCoefficientLog(srLength, i)+i*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE;
			}
//			System.out.println(Arrays.toString(logBinomD));
			logBinomialDesnity.put(srLength, logBinomD);
		}
//		System.out.println(maxDist);

//		int dist = LikelihoodUtils.hamDist(reads, subH);
//		double logProb = -plambda + dist * logPlambda - allFactorialLog[dist];
//		double logProb = ArithmeticUtils.binomialCoefficientLog(srLength, dist)+dist*LOG_ERROR_RATE+(srLength-dist)*LOG_ONE_MINUS_ERROR_RATE;

	}


    
	@Override
	public double getLogLikelihood(){
//		likelihoodKnown = false; //TODO REMOVE later //
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
        
	    
//	    double logLikelihood = calculateShoreReadLikelihoodBinomialModel2();
	    
//	    timeTrial();
//		return logLikelihood;
	}
	
	protected double calculateLogLikelihood() {
		double logLikelihood = calculateLogLikelihoodSelect(0);
		return logLikelihood;
	}

	public double calculateLogLikelihoodSelect(int index){

		double logLikelihood = Double.NEGATIVE_INFINITY;
		
		if (!Arrays.equals(NULL_SWAPINFO,  swapInfo) ){
			index=1;
		}
		switch (index) {
		case 0:
			
			logLikelihood = calculateShoreReadLikelihoodBinomalFull();
			break;
		case 1:
			
			logLikelihood = calculateShoreReadLikelihoodBinomal2();
			
			break;
//		case 2:
//			logLikelihood = calculateShoreReadLikelihoodBinomialModel2();
//			break;
		default:
			break;
		}
//	    double logLikelihood = calculateShoreReadLikelihood4();
//	    double logLikelihood = calculateShoreReadLikelihoodBinomialModel2();
	    
//	    timeTrial();
//		storeState();
		return logLikelihood;
	}


	private double calculateShoreReadLikelihoodBinomalFull() {
//System.out.println("FULL Likelihood calculation");
		double logLikelihood = 0;
		LikelihoodScaler liS = new LikelihoodScaler(LOG_C);
		
		int[] counter = new int[maxDist];
		for (int i = 0; i < aMap.getSrpCount(); i++) {

			Arrays.fill(counter, 0);
			char[] srpChar = shortReadChars.get(i);
			int start = aMap.getSrpStart(i);
			int end = aMap.getSrpEnd(i);
//	System.out.println("=========\ni_start_end\t"+i +"\t"+ start +"\t"+ end);
//			for (char[] hapCharArray : haplotypesChars) {
			for (int j = 0; j < haplotypesChars.length; j++) {
				
//				char[] hapCharArray = haplotypesChars[j]; 

				int dists = LikelihoodUtils.Dist(start, end, srpChar, haplotypesChars[j]);
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
			eachLikelihood[i] = liS.getLogLikelihood();
//			tempEachLikelihood[i] =eachLikelihood[i] ;
			logLikelihood += eachLikelihood[i];
			liS.reset(LOG_C);
		}	
//		storeState();
		
		
		return logLikelihood;
	}

	private double calculateShoreReadLikelihoodBinomal2() {
//		System.out.println("Partial");
//		System.out
//		.println("\n==========\nEQUAL Likelihood:" +"\t"+ Arrays.equals(eachLikelihood, tempEachLikelihood));
//		for (int i = 0; i < tempEachLikelihood.length; i++) {
//			System.out.println(i +"\t"+ tempEachLikelihood[i] +"\t"+ eachLikelihood[i]);
//		}
		
		double logLikelihood = 0;
		LikelihoodScaler liS = new LikelihoodScaler(LOG_C);
		
		int swapPos = swapInfo[1];
		ArrayList<Integer> mapPos = aMap.getMapToSrp(swapPos);
		
		
		int[] counter = new int[maxDist];
//		
//		double[] tempL = new double[aMap.getSrpCount()];
//		for (int i = 0; i < aMap.getSrpCount(); i++) {
//			Arrays.fill(counter, 0);
//			char[] srpChar = shortReadChars.get(i);
//			int start = aMap.getStart(i);
//			int end = aMap.getEnd(i);
//
//			for (int j = 0; j < haplotypesChars.length; j++) {
//				int dists = LikelihoodUtils.Dist(start, end, srpChar, haplotypesChars[j]);
//				counter[dists]++;
//			}	
//			
//			double[] logPD = logBinomialDesnity.get(aMap.getSrpLength(i));
//			for (int j = 0; j < counter.length; j++) {
//				if(counter[j]!=0){
//					liS.scaleLogProb(logPD[j], counter[j]);
//				}
//			}
//
//			tempL[i] = liS.getLogLikelihood();
////			if (mapPos.contains(i)){
////				System.out.println("add\t"+i);
////				tempEachLikelihood[i] = liS.getLogLikelihood();
////			}
//			liS.reset(LOG_C);
//		}	
		
//		for (int i = 0; i < aMap.getSrpCount(); i++) {
		for (int i : mapPos) {
					
			Arrays.fill(counter, 0);
			char[] srpChar = shortReadChars.get(i);
			int start = aMap.getSrpStart(i);
			int end = aMap.getSrpEnd(i);

			for (int j = 0; j < haplotypesChars.length; j++) {
				int dists = LikelihoodUtils.Dist(start, end, srpChar, haplotypesChars[j]);
				counter[dists]++;
			}	
			
			double[] logPD = logBinomialDesnity.get(aMap.getSrpLength(i));
			for (int j = 0; j < counter.length; j++) {
				if(counter[j]!=0){
					liS.scaleLogProb(logPD[j], counter[j]);
				}
			}

			eachLikelihood[i] = liS.getLogLikelihood();
			liS.reset(LOG_C);
		}

		logLikelihood = StatUtils.sum(eachLikelihood);
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

	
	
	
	@Override
	public Model getModel() {
		return this;
		
	}



	@Override
	public void makeDirty() {
        likelihoodKnown = false;
		
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == alignmentModel) {
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

	@Override
	protected void storeState() {
//		System.err.println("SR likelihood store: " + logLikelihood);
		System.arraycopy(eachLikelihood, 0, storedEachLikelihood, 0, eachLikelihood.length);
		storedLogLikelihood = logLikelihood;
//		System.arraycopy(storedEachLikelihood, 0, eachLikelihood, 0, eachLikelihood.length);
		
	}

	@Override
	public void restoreState() {
//		System.err.println("SR likelihood restore: " + storedLogLikelihood);
		logLikelihood = storedLogLikelihood;
		System.arraycopy(storedEachLikelihood, 0, eachLikelihood, 0, eachLikelihood.length);
		
	}

	@Override
	public void acceptState() {
//		System.out.println("Accept");
		
//		System.arraycopy(storedEachLikelihood, 0, eachLikelihood, 0, eachLikelihood.length);
		
	}

	@Deprecated
	public ShortReadLikelihood(Sequences reads, Alignment alignment){
		this(SHORT_READ_LIKELIHOOD);
	
		ArrayList<String> haplotypes = new ArrayList<String>();
		ArrayList<String> shortRead = new ArrayList<String>();
		for (int i = 0; i < alignment.getSequenceCount(); i++) {
			haplotypes.add(alignment.getAlignedSequenceString(i));
		}
		for (int i = 0; i < reads.getSequenceCount(); i++) {
			shortRead.add(reads.getSequence(i).getSequenceString());
		}
	
		setHaplotypes(haplotypes);
		setShortRead(shortRead);
	}

	@Deprecated
		public ShortReadLikelihood(ArrayList<String> shortRead, ArrayList<String> haplotypes){
			this(SHORT_READ_LIKELIHOOD);
			setHaplotypes(haplotypes);
			
			setShortRead(shortRead);
			
	//		preprocessHaplotypes();
	//		preprocessLikelihood();
		}

	/**
	 * @param shortRead the shortRead to set
	 */
	@Deprecated
	public void setShortRead(ArrayList<String> shortRead) {
		this.shortRead = shortRead;
		preprocessLikelihood();
	}

	/**
	 * @param haplotypes the haplotypes to set
	 */
	@Deprecated
	public void setHaplotypes(ArrayList<String> haplotypes) {
		this.haplotypes = haplotypes;
		this.hapLength = this.haplotypes.get(0).length();
		preprocessHaplotypes();
	}

	//	private void preprocessHaplotypes(char[][] charArray) {
	//		haplotypesCharsList.clear();
	//		for (String s : haplotypes) {
	//
	//			haplotypesCharsList.add(s.toCharArray());
	//		}
	//		
	//	}
		@Deprecated
		private void preprocessHaplotypes() {
			haplotypesCharsList.clear();
			for (String s : haplotypes) {
	
				haplotypesCharsList.add(s.toCharArray());
			}
			
		}

	//	public ShortReadLikelihood(ArrayList<String> shortRead, ArrayList<String> haplotypes,
	//			Parameter value) {
	//		this(shortRead, haplotypes);
	//		this.tempValue = value;
	//	}
	//
	//	public ShortReadLikelihood(ArrayList<String> shortRead, AlignmentModel alignmentModel,
	//			Parameter value) {
	//		//TODO FIXME
	//		this(SHORT_READ_LIKELIHOOD);
	//		setShortRead(shortRead);
	//        this.alignmentModel = alignmentModel;
	//		addModel(alignmentModel);
	//		this.tempValue = value;
	//	}
	//	
		@Deprecated
		private void preprocessLikelihood() {//~4ms
			
			shortReadChars = new ArrayList<char[]>();
	
			allFactorialLog = new double[hapLength+1];
			for (int i = 0; i < allFactorialLog.length; i++) {
				allFactorialLog[i] = ArithmeticUtils.factorialLog(i);
			}
			
			maxDist=0;
	
			logBinomialDesnity = new HashMap<Integer, double[]>();
			for (String reads : shortRead) {
				shortReadChars.add(reads.toCharArray());
				int srLength = reads.length();
	//			char[] srCharArray = reads.toCharArray();
				double plambda = srLength * ERROR_RATE;
				double logPlambda = Math.log(plambda);
	
				int srLength1 = srLength+1;
				maxDist = Math.max(maxDist, srLength1);
				double[] logBinomD = new double[srLength1];
				for (int i = 0; i < logBinomD.length; i++) {
	
					logBinomD[i] = ArithmeticUtils.binomialCoefficientLog(srLength, i)+i*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE;
				}
	
				logBinomialDesnity.put(srLength, logBinomD);
			}
			System.out.println(maxDist);
	
	//		int dist = LikelihoodUtils.hamDist(reads, subH);
	//		double logProb = -plambda + dist * logPlambda - allFactorialLog[dist];
	//		double logProb = ArithmeticUtils.binomialCoefficientLog(srLength, dist)+dist*LOG_ERROR_RATE+(srLength-dist)*LOG_ONE_MINUS_ERROR_RATE;
	
		}

	@Deprecated
	private ArrayList<String> shortRead;// = new ArrayList<>();
	@Deprecated
	private ArrayList<String> haplotypes;// = new ArrayList<>();
	@Deprecated
	private ArrayList<char[]> haplotypesCharsList = new ArrayList<char[]>();
	
}