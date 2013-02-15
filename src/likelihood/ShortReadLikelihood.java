package likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import model.AlignmentModel;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.util.ArithmeticUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import alignment.AlignmentMapping;
import alignment.AlignmentMatrix;

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
	
//	public static final int POISSON = 0;
	public static final int BINOMIAL = 0;
	
	private double logLikelihood;
    private double storedLogLikelihood;
    protected boolean likelihoodKnown = false;
	
	
	
	private ArrayList<String> shortRead;// = new ArrayList<>();
	private ArrayList<String> haplotypes;// = new ArrayList<>();
	private int maxDist;
	private double[] allFactorialLog;

	private ArrayList<char[]> haplotypesCharsList = new ArrayList<>();
	private char[][] haplotypesChars;
	private ArrayList<char[]> shortReadChars;
	
	
//	private HashMap<Integer, double[]> logPoissonDesnity;
	private HashMap<Integer, double[]> logBinomialDesnity;
	private Parameter tempValue;
	private AlignmentModel alignmentModel;
	private int hapLength;


	
	private AlignmentMapping aMap;
	
	
	
    @Override
	public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

    public ShortReadLikelihood(String name) {

        super(SHORT_READ_LIKELIHOOD);
        likelihoodKnown = false;
//
//        this.trialsParameter = trialsParameter;
//        this.proportionParameter = proportionParameter;
//        addVariable(trialsParameter);
//        addVariable(proportionParameter);
//        this.counts = counts;

    }

	public ShortReadLikelihood(AlignmentMapping aMap, AlignmentMatrix aMatrix){
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
		updatehaplotypesChars(aMatrix);
//		this.haplotypesChars = alignment.getCharMatrix();
		preprocessLikelihoodAlignmentMap();
	}
	public void updatehaplotypesChars(AlignmentMatrix aMatrix){
		this.haplotypesChars = aMatrix.getCharMatrix();
		makeDirty();
	}
	

	private void preprocessLikelihoodAlignmentMap() {//~4ms
		makeDirty();
		shortReadChars = new ArrayList<>();

		allFactorialLog = new double[hapLength+1];
		for (int i = 0; i < allFactorialLog.length; i++) {
			allFactorialLog[i] = ArithmeticUtils.factorialLog(i);
		}
		
		maxDist=0;

		logBinomialDesnity = new HashMap<>();
//		for (String reads : shortRead) {
		for (int s = 0; s < aMap.getSrpCount(); s++) {
			String srp = aMap.getFragmentSrp(s);
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


    
//	public ShortReadLikelihood(Sequences reads, SimpleAlignment alignment){
    public ShortReadLikelihood(Sequences reads, Alignment alignment){
		this(SHORT_READ_LIKELIHOOD);

		ArrayList<String> haplotypes = new ArrayList<>();
		ArrayList<String> shortRead = new ArrayList<>();
		for (int i = 0; i < alignment.getSequenceCount(); i++) {
			haplotypes.add(alignment.getAlignedSequenceString(i));
		}
		for (int i = 0; i < reads.getSequenceCount(); i++) {
			shortRead.add(reads.getSequence(i).getSequenceString());
		}

		setHaplotypes(haplotypes);
		setShortRead(shortRead);
	}


	private void preprocessHaplotypes(char[][] charArray) {
		haplotypesCharsList.clear();
		for (String s : haplotypes) {

			haplotypesCharsList.add(s.toCharArray());
		}
		
	}
	private void preprocessHaplotypes() {
		haplotypesCharsList.clear();
		for (String s : haplotypes) {

			haplotypesCharsList.add(s.toCharArray());
		}
		
	}

	public ShortReadLikelihood(ArrayList<String> shortRead, ArrayList<String> haplotypes){
		this(SHORT_READ_LIKELIHOOD);
		setHaplotypes(haplotypes);
		
		setShortRead(shortRead);
		
//		preprocessHaplotypes();
//		preprocessLikelihood();
	}

	public ShortReadLikelihood(ArrayList<String> shortRead, ArrayList<String> haplotypes,
			Parameter value) {
		this(shortRead, haplotypes);
		this.tempValue = value;
	}

	public ShortReadLikelihood(ArrayList<String> shortRead, AlignmentModel alignmentModel,
			Parameter value) {
		//TODO FIXME
		this(SHORT_READ_LIKELIHOOD);
		setShortRead(shortRead);
        this.alignmentModel = alignmentModel;
		addModel(alignmentModel);
		this.tempValue = value;
	}
	
	private void preprocessLikelihood() {//~4ms
		
		shortReadChars = new ArrayList<>();

		allFactorialLog = new double[hapLength+1];
		for (int i = 0; i < allFactorialLog.length; i++) {
			allFactorialLog[i] = ArithmeticUtils.factorialLog(i);
		}
		
		maxDist=0;

		logBinomialDesnity = new HashMap<>();
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

	@Override
	public double getLogLikelihood(){

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
		double logLikelihood = getLogLikelihoodSelect(0);
		return logLikelihood;
	}

	public double getLogLikelihoodSelect(int index){

		double logLikelihood = Double.NEGATIVE_INFINITY;
		switch (index) {
		case 0:
			logLikelihood = calculateShoreReadLikelihoodBinomal1();
			break;
//		case 1:
//			logLikelihood = calculateShoreReadLikelihood2();
//			break;
//		case 2:
//			logLikelihood = calculateShoreReadLikelihoodBinomialModel2();
//			break;
		default:
			break;
		}
//	    double logLikelihood = calculateShoreReadLikelihood4();
//	    double logLikelihood = calculateShoreReadLikelihoodBinomialModel2();
	    
//	    timeTrial();
		return logLikelihood;
	}


	private double calculateShoreReadLikelihoodBinomal1() {

		double logLikelihood = 0;
		LikelihoodScaler liS = new LikelihoodScaler(LOG_C);
		
		int[] counter = new int[maxDist];
	
		for (int i = 0; i < aMap.getSrpCount(); i++) {
			
		
			Arrays.fill(counter, 0);

			int start = aMap.getStart(i);
			int end = aMap.getEnd(i);
//	System.out.println("=========\ni_start_end\t"+i +"\t"+ start +"\t"+ end);
//			for (char[] hapCharArray : haplotypesChars) {
			for (int j = 0; j < haplotypesChars.length; j++) {
				
				char[] hapCharArray = haplotypesChars[j]; 

				int dists = LikelihoodUtils.Dist(start, end, shortReadChars.get(i), hapCharArray);
//				int noComb = hapCharArray.length - srCharArray.length + 1;
//				for (int j = 0; j < noComb; j++) {
				counter[dists]++;
//
//				}
			}	
			
			double[] logPD = logBinomialDesnity.get(aMap.getSrpLength(i));
			for (int j = 0; j < counter.length; j++) {
				if(counter[j]!=0){
//	System.out.println(j +" outOf "+(end-start) +"\t"+ logPD[j] +"\t"+ counter[j]);
					liS.scaleLogProb(logPD[j], counter[j]);
				}
			}
//	System.out.println("===============");
			logLikelihood += liS.getLogLikelihood();
			liS.reset(LOG_C);
		}	
		
		return logLikelihood;
	}

	public double calculateShoreReadLikelihoodBinomialModel2(){

		double logLikelihood = 0;
		LikelihoodScaler liS = new LikelihoodScaler(LOG_C);
		
		int[] counter = new int[maxDist];
	
		for (char[] srCharArray : shortReadChars){
			Arrays.fill(counter, 0);
						
			for (char[] hapCharArray : haplotypesCharsList) {

				int[] dists = LikelihoodUtils.hamDistAll(srCharArray, hapCharArray);
				int noComb = hapCharArray.length - srCharArray.length + 1;
				for (int j = 0; j < noComb; j++) {
					counter[dists[j]]++;

				}
			}	
			
			double[] logPD = logBinomialDesnity.get(srCharArray.length);
			for (int i = 0; i < counter.length; i++) {
				if(counter[i]!=0){
//					double logProb = -plambda +  logPD[i] - allFactorialLog[i];
					liS.scaleLogProb(logPD[i], counter[i]);
				}
			}

			logLikelihood += liS.getLogLikelihood();
			liS.reset(LOG_C);
		}	
		
		
		return logLikelihood;
		
	}
	
	public double calculateShoreReadLikelihoodBinomialModel(){

		double logLikelihood = 0;
		LikelihoodScaler liS = new LikelihoodScaler(LOG_C);
		for (String reads : shortRead) {
			int srLength = reads.length();
			double plambda = srLength * ERROR_RATE;
			double logPlambda = Math.log(plambda);

			for (String h : haplotypes) {
				int hLength = h.length();
				int noComb = hLength - srLength + 1;

				for (int j = 0; j < noComb; j++) {
					
					String subH = h.substring(j, j + srLength);
					int dist = LikelihoodUtils.hamDist(reads, subH);
//					double logProb = -plambda + dist * logPlambda - allFactorialLog[dist];
					double logProb = ArithmeticUtils.binomialCoefficientLog(srLength, dist)+dist*LOG_ERROR_RATE+(srLength-dist)*LOG_ONE_MINUS_ERROR_RATE;
					liS.scaleLogProb(logProb);
				}
			
			}	
			logLikelihood += liS.getLogLikelihood();
			liS.reset(LOG_C);
		}
		return logLikelihood;
		
	}
	
	
	
	
	/**
	 * @param shortRead the shortRead to set
	 */
	public void setShortRead(ArrayList<String> shortRead) {
		this.shortRead = shortRead;
		preprocessLikelihood();
	}

	/**
	 * @param haplotypes the haplotypes to set
	 */
	public void setHaplotypes(ArrayList<String> haplotypes) {
		this.haplotypes = haplotypes;
		this.hapLength = this.haplotypes.get(0).length();
		preprocessHaplotypes();
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
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void restoreState() {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void acceptState() {
		// TODO Auto-generated method stub
		
	}
	
}