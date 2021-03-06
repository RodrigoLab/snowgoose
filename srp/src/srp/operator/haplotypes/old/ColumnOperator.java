package srp.operator.haplotypes.old;

import srp.evolution.haplotypes.old.OldHapOperation;
import srp.evolution.haplotypes.old.OldHaplotypeModel;
import srp.evolution.shortreads.AlignmentMapping;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;


public class ColumnOperator extends AbstractCoercableOperator {

	
	public final static String OPERATOR_NAME = ColumnOperator.class.getSimpleName();
	public final static OldHapOperation OP = OldHapOperation.SWAPCOLUMN;
	private OldHaplotypeModel haplotypeModel;
	private int noOfHap;
	private int haplotypeLength;
	private AlignmentMapping alignmentMapping;
	private int haplotypeCount;
	private Parameter frequency;

	
	
//	public AlignmentSwapBaseOperator(Parameter parameter, HaplotypeModel haplotypeModel, int index, CoercionMode mode) {
////		super(mode);
//	}

	
	public ColumnOperator(OldHaplotypeModel haplotypeModel, int noOfHap, Parameter freqs, CoercionMode mode) {
		super(mode);

		this.frequency = freqs;
		this.haplotypeModel = haplotypeModel;
		this.noOfHap = noOfHap;
		
		haplotypeLength = this.haplotypeModel.getHaplotypeLength();
		haplotypeCount = this.haplotypeModel.getHaplotypeCount();
		this.noOfHap = haplotypeCount;
				
		alignmentMapping = this.haplotypeModel.getAlignmentMapping();
		
		
//		allPosChars = new int[2][haplotypeLength];
	
		
//		checkParameterIsValid();
		
//		scaleFactor = (int) (haplotypeLength*0.01);

//		if (scaleFactor <1) {
//			scaleFactor = 1;
//		}
		
//		convertToAutoOptimize(this.noOfHap);
	}

	@Override
	public String getPerformanceSuggestion() {

		return "";
	}

	@Override
	public String getOperatorName() {

		return OPERATOR_NAME;
	}


    @Override
	public double doOperation() throws OperatorFailedException {
		

		haplotypeModel.startHaplotypeOperation();

//		int[] posChar = alignmentMapping.getNextBaseFrequency(frequency);
		int[] posChar = haplotypeModel.getNextBaseFrequency(frequency);
		
		haplotypeModel.swapHaplotypeColumn(posChar);
		
		haplotypeModel.endHaplotypeOperation();
		
		return 0.0;//Incorrect!!??
		
		
		
	}

	@Override
	public double getCoercableParameter() {

		return 0;
	}

	@Override
	public void setCoercableParameter(double value) {
		
	}

	@Override
	public double getRawParameter() {
		return 0;
	}

}
