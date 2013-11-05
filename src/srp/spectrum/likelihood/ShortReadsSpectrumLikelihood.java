package srp.spectrum.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.ArithmeticUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.Haplotype;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.Operation;
import srp.haplotypes.ShortRead;
import srp.haplotypes.SwapInfo;
import srp.likelihood.LikelihoodScaler;
import srp.likelihood.LikelihoodUtils;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import dr.evolution.datatype.DataType;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

public class ShortReadsSpectrumLikelihood  extends AbstractModelLikelihood {

	
    public static final String SHORT_READ_LIKELIHOOD = "ShortReadLikelihood";
	public static final String NAME = SHORT_READ_LIKELIHOOD;
	public static final double ERROR_RATE = 0.0107;
	public static final double NOT_ERROR_RATE = 1-ERROR_RATE;
	public static final double LOG_ERROR_RATE = Math.log(ERROR_RATE);
	public static final double LOG_ONE_MINUS_ERROR_RATE = Math.log(1-ERROR_RATE);
	public static final double C = 1e-200;
	public static final double LOG_C = Math.log(C);
	
	private static final double EVALUATION_TEST_THRESHOLD = 1e-8;
	//	public static final int[] NULL_SWAPINFO = HaplotypeModel.NULL_SWAPINFO;
	protected boolean likelihoodKnown = false;
	
	private int spectrumLength;
	private int spectrumCount;
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
	
//	private HaplotypeModel haplotypeModel;
	private LikelihoodScaler liS;
	private SwapInfo swapInfo;
	
	@Deprecated
	private int[] counter;
	private Operation operation;
	private SpectrumAlignmentModel spectrumModel;
	private DataType dataType;
	
	
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
		this.swapInfo = this.spectrumModel.getSwapInfo();
		this.dataType = this.spectrumModel.getDataType();
		operation = Operation.NONE;
		preprocessLikelihoodAlignmentMap();
		calculateSrpLikelihoodFull();
//		calculateSrpLikelihoodFullUseCounter();
		
		addModel(this.spectrumModel);
	}
	

	private void preprocessLikelihoodAlignmentMap() {
		makeDirty();
		
		liS = new LikelihoodScaler(LOG_C);
		
		srpCount = aMap.getSrpCount();
		spectrumCount = spectrumModel.getSpectrumCount();
		spectrumLength = spectrumModel.getSpectrumLength();
		
		logLikelihood = Double.NEGATIVE_INFINITY;
		storedLogLikelihood = Double.NEGATIVE_INFINITY;
		
		allDists = new int[srpCount][spectrumCount];
		storedAllDists = new int[srpCount][spectrumCount];

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
		SwapInfo swapInfo = spectrumModel.getSwapInfo();
		operation = swapInfo.getOperation();

		double logLikelihood = Double.NEGATIVE_INFINITY;

		switch (operation) {
			case NONE:
				logLikelihood = calculateSrpLikelihoodFull();
				break;
//			case SWAPMULTI:
//				logLikelihood = calculateSrpLikelihoodMultiBasesSwap();
//				break;

			case SWAPSINGLE:
				logLikelihood = calculateSrpLikelihoodSingleBaseSwap(swapInfo);
				break;
//			case UNIFORMSWAPBASE:
//				logLikelihood = calculateSrpLikelihoodSingleBaseSwap();
//				break;
//			case SWAPCOLUMN:
//				logLikelihood = calculateSrpLikelihoodSwapColumn();
//				break;
//			case SWAPSECTION:
//				logLikelihood = calculateSrpLikelihoodSwapSection();
//				break;
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
			String fullSrp = aMap.getSrpFull(i);
			int start = aMap.getSrpStart(i);
			int end = aMap.getSrpEnd(i);
			
			double[] logPD = scaledLogBinomialDesnity.get(aMap.getSrpLength(i));

			liS.reset();
			for (int j = 0; j < spectrumCount; j++) {

				Spectrum spectrum = spectrumModel.getSpectrum(j);
				double spectrumLogLikelihood = 0;
				for (int k = start; k < end; k++) {
					double[] frequencies = spectrum.getFrequencies(k);
					char srpChar = fullSrp.charAt(k);
					int state = dataType.getState(srpChar);
					double likelihood = frequencies[state] * NOT_ERROR_RATE
							+ (1 - frequencies[state]) * ERROR_RATE;
					double logLikelihood = Math.log(likelihood);
					System.out.println(likelihood +"\t"+ Arrays.toString(frequencies));
//					for (int l = 0; l < frequencies.length; l++) {
//						
//					}
					spectrumLogLikelihood += logLikelihood;
					
				}
				System.out.println("SL: "+spectrumLogLikelihood);
				liS.scaleLogProb(spectrumLogLikelihood);
				
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
		
		return logLikelihood;
	}


	private double calculateSrpLikelihoodSingleBaseSwap(SwapInfo swapInfo) {
		
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
			for (int j = 0; j < spectrumCount ; j++) {
				liS.addScaledLogProb(logPD[allDists[srpIndex][j]]);
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
        if (model == spectrumModel) {
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
	}

	public Operation getOperation(){
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

}