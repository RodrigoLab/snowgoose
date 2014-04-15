package srp.likelihood;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.carrotsearch.hppc.BitSet;

import srp.dr.evolution.datatype.ShortReads;
import srp.evolution.AbstractAlignmentModel;
import srp.evolution.OperationRecord;
import srp.evolution.OperationType;
import srp.evolution.shortreads.ShortReadMapping;
import srp.evolution.spectrum.AbstractSpectrumAlignmentModel;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import srp.haplotypes.AbstractHaplotypeModel;
import srp.likelihood.stateLikelihood.BetaMeanStateLikelihood;
import srp.likelihood.stateLikelihood.BetaModeStateLikelihood;
import srp.likelihood.stateLikelihood.ChisqStateLikelihood;
import srp.likelihood.stateLikelihood.GTestStateLikelihood;
import srp.likelihood.stateLikelihood.ProbabilityStateLikelihood;
import srp.likelihood.stateLikelihood.StateLikelihood;
import dr.evolution.datatype.DataType;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

public abstract class AbstractShortReadsLikelihood extends
		AbstractModelLikelihood {

	private static final long serialVersionUID = 2079474866153379297L;
	public static final double ERROR_RATE = 0.0107;
	public static final double NOT_ERROR_RATE = 1-ERROR_RATE;
	public static final double LOG_ERROR_RATE = Math.log(ERROR_RATE);
	public static final double LOG_NOT_ERROR_RATE = Math.log(NOT_ERROR_RATE);
	public static final double LOG_ONE_MINUS_ERROR_RATE = Math.log(1-ERROR_RATE);
	public static final double C = 1e-300;
	public static final double LOG_C = Math.log(C);

	public static final DataType DATA_TYPE = ShortReads.INSTANCE;
	public static final int STATE_COUNT = DATA_TYPE.getStateCount();//4
	public static final int AMBIGUOUS_STATE_COUNT = DATA_TYPE.getAmbiguousStateCount();//18
	public static final int GAP_STATE = DATA_TYPE.getGapState();;
	
	private static final double EVALUATION_TEST_THRESHOLD = 1e-8;
	public static final String SHORT_READ_LIKELIHOOD = "ShortReadsLikelihood";

	protected boolean debug = false;;
	
	protected double logLikelihood;
	protected double storedLogLikelihood;
	
	protected double[] eachSrpLogLikelihood;
	protected double[] storedEachSrpLogLikelihood;
	
	protected double[] sumScaledSrpLogLikelihood;
	protected double[] storedSumScaledSrpLogLikelihood;


	protected boolean likelihoodKnown;
	
	protected int sequenceLength;
	protected int sequenceCount;
	protected int srpCount;

	protected LikelihoodScaler liS;
	protected ShortReadMapping srpMap;
	protected OperationRecord operationRecord;

	protected MultiType multiType;
	protected boolean[] srpSwitch; //MultiType.Array
	protected Set<Integer> allSrpPos; //MultiType.Hash
	protected BitSet bitSet; //MultiType.BitSet
	protected int[] srpIndex; //MultiType.BitSet
	protected int srpIndexCount; //MultiType.BitSet

	protected int[][] mapToSrpArray;	
	
	protected int[][] allSrpState2D;
	protected char[][] allSrpChar2D;
//	protected String[] srpArray;

	public AbstractShortReadsLikelihood(String name, ShortReadMapping srpMap) {
		super(name);
		liS = new LikelihoodScaler(LOG_C);
		preprocessShortReadMapping(srpMap);
	}


	private void preprocessShortReadMapping(ShortReadMapping srpMap2) {

		this.srpMap = srpMap2;
		srpCount = srpMap.getSrpCount();
		
		srpSwitch = new boolean[srpCount];
		allSrpPos = new HashSet<Integer>();
		bitSet = new BitSet(srpCount);
		srpIndex = new int[srpCount];
		
		logLikelihood = Double.NEGATIVE_INFINITY;
		storedLogLikelihood = Double.NEGATIVE_INFINITY;

		sumScaledSrpLogLikelihood = new double[srpCount];
		storedSumScaledSrpLogLikelihood = new double[srpCount];
		
		eachSrpLogLikelihood = new double[srpCount];
		storedEachSrpLogLikelihood = new double[srpCount];

		
		mapToSrpArray = srpMap.getMapToSrpArray();
		allSrpState2D = srpMap.getSrpState2DArray();
		allSrpChar2D = srpMap.getSrpChar2DArray();
				
		
		
	}


	protected double calculateLogLikelihood() {

//		OperationType operation = alignmentModel.getOperation();
		OperationType operation = operationRecord.getOperation();
		double logLikelihood = Double.NEGATIVE_INFINITY;

		if(debug){
			System.out.println("Calculate ShortReadLikelihood:\t"+operation);
		}
		switch (operation) {
		case NONE:
//			logLikelihood = calculateSrpLikelihoodFull();
			logLikelihood = calculateSrpLikelihoodFull();
			break;
		case FULL:
			logLikelihood = calculateSrpLikelihoodFull();
			break;
		case SINGLE:
			logLikelihood = calculateSrpLikelihoodSingle();
			break;
		case COLUMN:
			logLikelihood = calculateSrpLikelihoodColumn();
			break;
		case MULTI:
			logLikelihood = calculateSrpLikelihoodMulti();
			break;
		case SWAP_SUBCOLUMN:
			logLikelihood = calculateSrpLikelihoodSubColumn();
			break;
		case RECOMBINATION:
			logLikelihood = calculateSrpLikelihoodRecombination();
			break;
		// case MASTER:
		// logLikelihood = calculateSrpLikelihoodFullMaster()
		// break;
		default:
			throw new IllegalArgumentException("Unknown operation type: "
					+ operation);

		}

		return logLikelihood;
	}
	
	protected abstract double calculateSrpLikelihoodRecombination();

	protected abstract double calculateSrpLikelihoodSubColumn();

	protected abstract double calculateSrpLikelihoodMulti();

	protected abstract double calculateSrpLikelihoodColumn();

	protected abstract double calculateSrpLikelihoodSingle();

	protected abstract double calculateSrpLikelihoodFullMaster();

	protected abstract double calculateSrpLikelihoodFull();


	protected double getStoredLogLikelihood() {
		return storedLogLikelihood;
	}
	


	protected void recalculateArray(int[] siteIndexs) {
		for (int s : siteIndexs) {
			for (int i : mapToSrpArray[s]){
				srpSwitch[i] = true;
			}
		}
	}

	protected void recalculateHashSet(int[] siteIndexs) {
		allSrpPos.clear();
		for (int s : siteIndexs) {
			ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(s);
			allSrpPos.addAll(mapToSrp);
		}
	}

	protected void recalculateBitSet(int[] siteIndexs) {
		bitSet.clear();
//		BitSet bitSet = new BitSet(srpCount);
		for (int s : siteIndexs) {
			BitSet tempSet = srpMap.getBitSet(s);
			bitSet.or(tempSet);
		}
	}


	public void setMultiType(MultiType type) {
		this.multiType = type;
	}

	public OperationType getOperation() {
		return operationRecord.getOperation();
	}
	
	protected static int getStateAtK(String fullSrp, int k) {
		char srpChar = fullSrp.charAt(k);//TODO: change to char[] at hight lv or at ShortReadMapping
		int state = DATA_TYPE.getState(srpChar);
		
		return state;
	}


	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
        
        likelihoodKnown = false;
//		if(model instanceof AbstractHaplotypeModel){
//        if(model instanceof AbstractAlignmentModel){
//			operationRecord = ((AbstractAlignmentModel) model).getOperationRecord();
//			System.out.println(operationRecord.getOperation());
//		}
	}

	@SuppressWarnings("rawtypes")
	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		System.err.println("Call handleVariableChangedEvent in SpectrumAlignmentModel");
	}

	@Override
	public double getLogLikelihood() {
	
		if (!likelihoodKnown) {
			logLikelihood = calculateLogLikelihood();
			likelihoodKnown = true;
		}
		return logLikelihood;
	}


	@Override
	public Model getModel() {
		return this;
	
	}

	@Override
	public Element createElement(Document d) {
		throw new RuntimeException("Not implemented yet!");
	}

	@Override
	protected void acceptState() {
		//Do nothing
	}

	protected enum MultiType {
		Array, Hash, All, BitSet, ;

	}

}