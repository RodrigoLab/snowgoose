package srp.likelihood;

import java.util.ArrayList;
import java.util.Set;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.carrotsearch.hppc.BitSet;

import srp.dr.evolution.datatype.ShortReads;
import srp.likelihood.stateLikelihood.BetaMeanStateLikelihood;
import srp.likelihood.stateLikelihood.BetaModeStateLikelihood;
import srp.likelihood.stateLikelihood.ChisqStateLikelihood;
import srp.likelihood.stateLikelihood.GTestStateLikelihood;
import srp.likelihood.stateLikelihood.ProbabilityStateLikelihood;
import srp.likelihood.stateLikelihood.StateLikelihood;
import srp.shortreads.ShortReadMapping;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
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
//	public static final double LOG_ONE_MINUS_ERROR_RATE = Math.log(1-ERROR_RATE);
	public static final double C = 1e-250;
	public static final double LOG_C = Math.log(C);

	public static final DataType DATA_TYPE = ShortReads.INSTANCE;
	public static final int STATE_COUNT = DATA_TYPE.getStateCount();//4
	public static final int AMBIGUOUS_STATE_COUNT = DATA_TYPE.getAmbiguousStateCount();//18
	
	private static final double EVALUATION_TEST_THRESHOLD = 1e-8;
	private static final int GAP_STATE = 17;  


	protected double logLikelihood;
	protected double storedLogLikelihood;
	
	protected boolean likelihoodKnown;
	
	protected int spectrumLength;
	protected int spectrumCount;

	protected MultiType multiType;

	protected ShortReadMapping srpMap;
	
	protected int[][] mapToSrpArray;
	
	protected boolean[] srpSwitch;
	protected Set<Integer> allSrpPos;
	protected BitSet bitSet;
	
	public AbstractShortReadsLikelihood(String name) {
		super(name);
	}

	protected abstract double calculateLogLikelihood();

	
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

	
	protected static int getStateAtK(String fullSrp, int k) {
		char srpChar = fullSrp.charAt(k);//TODO: change to char[] at hight lv or at ShortReadMapping
		int state = DATA_TYPE.getState(srpChar);
		
		return state;
	}


	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
        
        likelihoodKnown = false;
		
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