package likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import model.AlignmentModel;

import org.apache.commons.math3.util.ArithmeticUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import dr.evolution.alignment.Alignment;
import dr.evolution.sequence.Sequences;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;


@SuppressWarnings({ "unused", "serial" })
@Deprecated
public class ZOLD_ShortReadLikelihood extends AbstractModelLikelihood {

    public static final String SHORT_READ_LIKELIHOOD = "ShortReadLikelihood";
	public static final String NAME = SHORT_READ_LIKELIHOOD;
	public static final double ERROR_RATE = 0.0107;
	public static final double ONE_MINUS_ERROR_RATE = 1-ERROR_RATE;
	public static final double LOG_ERROR_RATE = Math.log(ERROR_RATE);
	public static final double LOG_ONE_MINUS_ERROR_RATE = Math.log(1-ERROR_RATE);
	public static final double C = 1e-200;
	public static final double LOG_C = Math.log(C);
	
	public static final int POISSON = 0;
	public static final int BINOMIAL = 2;
	
	private double logLikelihood;
    private double storedLogLikelihood;
    protected boolean likelihoodKnown = false;
	
	
	
	private ArrayList<String> shortRead;// = new ArrayList<>();
	private ArrayList<String> haplotypes;// = new ArrayList<>();
	private int maxDist;
	private double[] allFactorialLog;

	private ArrayList<char[]> haplotypesChars = new ArrayList<char[]>();
	private ArrayList<char[]> shortReadChars;
	
	private HashMap<Integer, double[]> logPoissonDesnity;
	private HashMap<Integer, double[]> logBinomialDesnity;
	private Parameter tempValue;
	private AlignmentModel alignmentModel;
	private int hapLength;
	
	
	
    @Override
	public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

    public ZOLD_ShortReadLikelihood(String name) {

        super(SHORT_READ_LIKELIHOOD);
        likelihoodKnown = false;
//
//        this.trialsParameter = trialsParameter;
//        this.proportionParameter = proportionParameter;
//        addVariable(trialsParameter);
//        addVariable(proportionParameter);
//        this.counts = counts;

    }

//	public ShortReadLikelihood(Sequences reads, SimpleAlignment alignment){
    public ZOLD_ShortReadLikelihood(Sequences reads, Alignment alignment){
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


	private void preprocessHaplotypes() {
		haplotypesChars.clear();
		for (String s : haplotypes) {
			haplotypesChars.add(s.toCharArray());
		}
		
	}

	public ZOLD_ShortReadLikelihood(ArrayList<String> shortRead, ArrayList<String> haplotypes){
		this(SHORT_READ_LIKELIHOOD);
		setHaplotypes(haplotypes);
		
		setShortRead(shortRead);
		
//		preprocessHaplotypes();
//		preprocessLikelihood();
	}

	public ZOLD_ShortReadLikelihood(ArrayList<String> shortRead, ArrayList<String> haplotypes,
			Parameter value) {
		this(shortRead, haplotypes);
		this.tempValue = value;
	}

	public ZOLD_ShortReadLikelihood(ArrayList<String> shortRead, AlignmentModel alignmentModel,
			Parameter value) {
		//TODO FIXME
		this(SHORT_READ_LIKELIHOOD);
		setShortRead(shortRead);
        this.alignmentModel = alignmentModel;
		addModel(alignmentModel);
		this.tempValue = value;
	}
	
	private void preprocessLikelihood() {//~4ms
		
		shortReadChars = new ArrayList<char[]>();

		allFactorialLog = new double[hapLength+1];
		for (int i = 0; i < allFactorialLog.length; i++) {
			allFactorialLog[i] = ArithmeticUtils.factorialLog(i);
		}
		
		maxDist=0;
		logPoissonDesnity = new HashMap<Integer, double[]>();
		logBinomialDesnity = new HashMap<Integer, double[]>();
		for (String reads : shortRead) {
			shortReadChars.add(reads.toCharArray());
			int srLength = reads.length();
//			char[] srCharArray = reads.toCharArray();
			double plambda = srLength * ERROR_RATE;
			double logPlambda = Math.log(plambda);

			int srLength1 = srLength+1;
			maxDist = Math.max(maxDist, srLength1);
			double[] logPoisD = new double[srLength1];
			double[] logBinomD = new double[srLength1];
			for (int i = 0; i < logPoisD.length; i++) {
//				logPD[i] = logPlambda*i;
				logPoisD[i] = -plambda + logPlambda*i -allFactorialLog[i];
				logBinomD[i] = ArithmeticUtils.binomialCoefficientLog(srLength, i)+i*LOG_ERROR_RATE+(srLength-i)*LOG_ONE_MINUS_ERROR_RATE;
			}
			logPoissonDesnity.put(srLength, logPoisD);
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
		double logLikelihood = calculateShoreReadLikelihood5();
		return logLikelihood;
	}

	public double getLogLikelihoodSelect(int index){

		double logLikelihood = Double.NEGATIVE_INFINITY;
		switch (index) {
		case 0:
			logLikelihood = calculateShoreReadLikelihood5();
			break;
		case 1:
			logLikelihood = calculateShoreReadLikelihood2();
			break;
		case 2:
			logLikelihood = calculateShoreReadLikelihoodBinomialModel2();
			break;
		default:
			break;
		}
//	    double logLikelihood = calculateShoreReadLikelihood4();
//	    double logLikelihood = calculateShoreReadLikelihoodBinomialModel2();
	    
//	    timeTrial();
		return logLikelihood;
	}


	private double calculateShoreReadLikelihood5(){

		double logLikelihood = 0;
		LikelihoodScaler liS = new LikelihoodScaler(LOG_C);
		
		int[] counter = new int[maxDist];
	
		for (char[] srCharArray : shortReadChars){
			Arrays.fill(counter, 0);
			
//			double plambda = srLength * ERROR_RATE;
//			double logPlambda = Math.log(plambda);
			
			for (String hapString : haplotypes) {
System.out.print(hapString.length() +"\t" );
				int[] dists = LikelihoodUtils.hamDistAll(srCharArray, hapString);
				int noComb = hapString.length() - srCharArray.length + 1;

				for (int j = 0; j < noComb; j++) {
					
					counter[dists[j]]++;
					
//					int d = dists[j];1
//					double logProb = -plambda +  logPD[d] - allFactorialLog[d];
//					double logProb = logPD[d];

//					liS.scaleLogProb(logProb);
	
				}
			}	
//			System.out.println(Arrays.toString(counter));
			double[] logPD = logPoissonDesnity.get(srCharArray.length);
			for (int i = 0; i < counter.length; i++) {
				if(counter[i]!=0){
//					double logProb = -plambda +  logPD[i] - allFactorialLog[i];
					liS.scaleLogProb(logPD[i], counter[i]);
					
				}
			}

			logLikelihood += liS.getLogLikelihood();
			liS.reset(LOG_C);
		}	
		System.out.println("Calculate shoreRdeaLikelihood out" +"\t"+ logLikelihood);
		return logLikelihood;
	}
	
	

	private double calculateShoreReadLikelihood4(){

		double logLikelihood = 0;
		LikelihoodScaler liS = new LikelihoodScaler(LOG_C);
		
		int[] counter = new int[maxDist];
	
		for (char[] srCharArray : shortReadChars){
			Arrays.fill(counter, 0);
			
//			double plambda = srLength * ERROR_RATE;
//			double logPlambda = Math.log(plambda);
			
			for (char[] hapCharArray : haplotypesChars) {

				int[] dists = LikelihoodUtils.hamDistAll(srCharArray, hapCharArray);
				int noComb = hapCharArray.length - srCharArray.length + 1;

				for (int j = 0; j < noComb; j++) {
					
					counter[dists[j]]++;
					
//					int d = dists[j];
//					double logProb = -plambda +  logPD[d] - allFactorialLog[d];
//					double logProb = logPD[d];

//					liS.scaleLogProb(logProb);
	
				}
			}	
			
			double[] logPD = logPoissonDesnity.get(srCharArray.length);
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
	
	
	private double calculateShoreReadLikelihood3(){

		double logLikelihood = 0;
		LikelihoodScaler liS = new LikelihoodScaler(LOG_C);
						
		for (String reads : shortRead) {
			int srLength = reads.length();
			char[] srCharArray = reads.toCharArray();
			double plambda = srLength * ERROR_RATE;
//			double logPlambda = Math.log(plambda);
			double[] logPD = logPoissonDesnity.get(srLength);

			for (char[] hapCharArray : haplotypesChars) {
				
				int[] dists = LikelihoodUtils.hamDistAll(srCharArray, hapCharArray);//, j+srLength);
				for (int j = 0; j < dists.length; j++) {
					int d = dists[j];
					
					double logProb = -plambda +  logPD[d] - allFactorialLog[d];
//					double logProb = -plambda +  d * logPlambda - allFactorialLog[d];

					liS.scaleLogProb(logProb);
	
				}
			}	
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
						
			for (char[] hapCharArray : haplotypesChars) {

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
	
	
	
	private double calculateShoreReadLikelihood2(){

		double logLikelihood = 0;
		LikelihoodScaler liS = new LikelihoodScaler(LOG_C);
		for (String reads : shortRead) {
			int srLength = reads.length();
			
			
			double plambda = srLength * ERROR_RATE;
			double logPlambda = Math.log(plambda);

			for (String h : haplotypes) {
				int hLength = h.length();
				int noComb = hLength - srLength + 1;
	int[] allDists = new int[noComb];
				for (int j = 0; j < noComb; j++) {
					
					String subH = h.substring(j, j + srLength);
					int dist = LikelihoodUtils.hamDist(reads, subH);
					allDists[j] = dist;
					double logProb = -plambda + dist * logPlambda - allFactorialLog[dist];
					liS.scaleLogProb(logProb);
				}

			}	
			logLikelihood += liS.getLogLikelihood();
			liS.reset(LOG_C);
		}
		return logLikelihood;
	}
	
	
	private double calculateShoreReadLikelihood1(){

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
					double logProb = -plambda + dist * logPlambda - ArithmeticUtils.factorialLog(dist);
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