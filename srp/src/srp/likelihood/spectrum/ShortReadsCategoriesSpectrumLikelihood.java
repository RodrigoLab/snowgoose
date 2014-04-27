package srp.likelihood.spectrum;

import java.util.BitSet;
import java.util.Set;

import srp.dr.evolution.datatype.ShortReads;
import srp.evolution.shortreads.ShortReadMapping;
import srp.evolution.spectrum.CategorySpectrumAlignmentModel;
import srp.evolution.spectrum.SpectraParameter;
import srp.likelihood.LikelihoodScaler;
import srp.likelihood.stateLikelihood.StateLikelihood;
import dr.evolution.datatype.DataType;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
//import java.util.BitSet;


public class ShortReadsCategoriesSpectrumLikelihood  extends AbstractShortReadsSpectrumLikelihood {

	private static final long serialVersionUID = 7438385718398999755L;

	private static final boolean DEBUG = false;
	
    public static final String SHORT_READ_LIKELIHOOD = "ShortReadSpectrumLikelihood";
	public static final String NAME = SHORT_READ_LIKELIHOOD;
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
	private static final int TWO = 2;
	private static final int GAP_STATE = 17;  
	

	
	private final double MIN_LOG_LIKELIHOOD;

	protected boolean likelihoodKnown;
	
	private int spectrumLength;
	private int spectrumCount;
	private int srpCount;

	private double logLikelihood;
	private double storedLogLikelihood;

	private double[] eachSrpLogLikelihood;
	private double[] storedEachSrpLogLikelihood;

	private double[] sumScaledSrpLogLikelihood;
	private double[] storedSumSrpLogLikelihood;
	
	
	private double[] spectrumLogLikelihood;
	private double[] storedSpectrumLogLikelihood;
	
	private double[] scaledSpectrumLogLikelihood;
	private double[] storedScaledSpectrumLogLikelihood;
	
	
	private double[] allStateLogLikelihood;
	private double[] allStoredStateLogLikelihood;
	
	private double[] allStateLogLikelihoodIJK;
	private double[] allStoredStateLogLikelihoodIJK;
	
	private StateLikelihood stateLikelihood;
	private ShortReadMapping srpMap;
	private CategorySpectrumAlignmentModel spectrumModel;
	private LikelihoodScaler liS;

	
	private int[] srpIndex ;
	private boolean[] srpSwitch;
	private Set<Integer> allSrpPos;

	private MultiType multiType;
//	private DistType distType;

	private String[] srpArray;
	private int[][] allSrpState2D;

	private int[][] mapToSrpArray;

	private BitSet bitSet;

	private int srpIndexCount;
	
//	public long timeStart = 0;//REMOVE
//	public long m1Time = 0;//REMOVE
//	public long m2Time = 0;//REMOVE
//	public long liTime  = 0;//REMOVE
//
//	double[][] m1State;// = new double[STATE_COUNT][spectrumLength];//REMOVE
//	double[][] m1StoreState;// = new double[STATE_COUNT][spectrumLength];//REMOVE
//	
//	double[][] m2State;// = new double[spectrumLength][STATE_COUNT];//REMOVE
//	double[][] m2StoreState;// = new double[spectrumLength][STATE_COUNT];//REMOVE
//	
//	double[] m3State;// = new double[STATE_COUNT][spectrumLength];//REMOVE
//	double[] m3StoreState;// = new double[STATE_COUNT][spectrumLength];//REMOVE

	
	
	public ShortReadsCategoriesSpectrumLikelihood(CategorySpectrumAlignmentModel spectrumModel, ShortReadMapping srpMap, DistType distType){
		
		super(SHORT_READ_LIKELIHOOD, srpMap);
		System.err.println("Warning!!! Incompletely implementation");
		this.spectrumModel = spectrumModel;


//		multiType = MultiType.Array;
		multiType = MultiType.BitSet;
		
//		type = MultiType.Hash;
//		type = MultiType.All;
//		distTypeCode = "flat";//"betaMean"  "betaMode" "gTest"
		setDistType(distType);
		MIN_LOG_LIKELIHOOD = stateLikelihood.caluclateStateLogLikelihood(SpectraParameter.MIN_FREQ);

		likelihoodKnown = false;
		System.err.println("Never finished implementing this class!!");
		System.exit(-1);
		
	}

	@Override
	protected double calculateLogLikelihood() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// TODO Auto-generated method stub
		
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

	@Override
	protected double calculateSrpLikelihoodRecombination() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected double calculateSrpLikelihoodSubColumn() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected double calculateSrpLikelihoodMulti() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected double calculateSrpLikelihoodColumn() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected double calculateSrpLikelihoodSingle() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected double calculateSrpLikelihoodFullMaster() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected double calculateSrpLikelihoodFull() {
		// TODO Auto-generated method stub
		return 0;
	}
	
}
