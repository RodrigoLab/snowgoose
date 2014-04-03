package srp.operator.haplotypes;

import java.util.Arrays;

import srp.haplotypes.old.OldHaplotypeModel;
import srp.haplotypes.old.OldHapOperation;
import srp.haplotypes.old.OldHapSwapInfo;
import srp.shortreads.AlignmentMapping;
import dr.evolution.datatype.Nucleotides;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class SwitchBaseFrequencyOperator extends AbstractCoercableOperator {

	public final static OldHapOperation OP = OldHapOperation.SWAPSINGLE;
	
//	protected final int haplotypeLength;

	protected OldHaplotypeModel haplotypeModel;

//	protected int swapLength;

//	protected int[][] allPosChars;
	protected AlignmentMapping alignmentMapping;

//	private double autoOptimize;
//	private double scaleFactor;

	public final static String OPERATOR_NAME = SwitchBaseFrequencyOperator.class.getSimpleName();
	
	private static final int NUCLEOTIDE_STATES[] = Nucleotides.NUCLEOTIDE_STATES;
	private static final int NEW_CHAR_INDEX = OldHapSwapInfo.SWAPBASE_NEW_CHAR_INDEX;
	private static final int OLD_CHAR_INDEX = OldHapSwapInfo.SWAPBASE_OLD_CHAR_INDEX;
	
	private Parameter frequency;

	private double switchProb;

	public SwitchBaseFrequencyOperator(OldHaplotypeModel haplotypeModel, double switchProb, Parameter freqs, CoercionMode mode )  {
		super(mode);
		this.frequency = freqs;
			
		this.haplotypeModel = haplotypeModel;
		this.switchProb = switchProb;
		
		setTargetAcceptanceProbability(targetAcceptanceProb);

		//		haplotypeLength = this.haplotypeModel.getHaplotypeLength();
//		alignmentMapping = this.haplotypeModel.getAlignmentMapping();
		
//		allPosChars = new int[2][haplotypeLength];
	
		
//		checkParameterIsValid();
		
//		scaleFactor = (int) (haplotypeLength*0.01);
//
//		if (scaleFactor <1) {
//			scaleFactor = 1;
//		}
//		
//		convertToAutoOptimize(this.swapLength);
		
	}

//    public double getCoercableParameter() {
//        return Math.log(delta);
//    }
//
//    public void setCoercableParameter(double value) {
//        delta = Math.exp(value);
//    }


	/*
	*
	* @return the log-transformed hastings ratio
	*/
	@Override
	public double doOperation() throws OperatorFailedException {

	
		double logq = 0.0;
		double rand = MathUtils.nextDouble();
		if(rand < switchProb){
			haplotypeModel.startHaplotypeOperation();
			logq = haplotypeModel.swapNextDiffBaseFrequency(OP, frequency);
			haplotypeModel.endHaplotypeOperation();
			count++;
		}
//		else{
//			haplotypeModel.storeOperationRecord(Operation.PASS);
//		}
		return logq;
	}
	int count = 0;
	@Override
	public double getCoercableParameter() {
	    return Math.log(switchProb / (1.0 - switchProb));
	}
	

	@Override
	public void setCoercableParameter(double value) {
		switchProb = Math.exp(value) / (1.0 + Math.exp(value));
	    
	}

	@Override
	public double getRawParameter() {
		return switchProb;
	}
	
	

	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

	@Override
	public String getPerformanceSuggestion() {
		return switchProb+ "\t"
				+ Arrays.toString(haplotypeModel.getSwapInfo()
						.getSwapInfoSWAPBASE()) + "\t"
				+ (double) getAcceptCount()
				/ (getRejectCount() + getAcceptCount()) + "::" + "\t"
				+ getCoercableParameter();

	}

    @Override
	public double getMinimumAcceptanceLevel() {
        return targetAcceptanceProb-0.2;
    }

    @Override
	public double getMaximumAcceptanceLevel() {
        return targetAcceptanceProb + 0.2;
    }

    @Override
	public double getMinimumGoodAcceptanceLevel() {
        return targetAcceptanceProb- 0.1;
    }

    @Override
	public double getMaximumGoodAcceptanceLevel() {
        return targetAcceptanceProb+ 0.1;
    }

    private double targetAcceptanceProb = 0.5;


	
}
