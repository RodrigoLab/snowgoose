package srp.operator.haplotypes;

import srp.evolution.OperationType;
import srp.haplotypes.Haplotype;
import srp.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;


public class ColumnOperator extends AbstractHaplotypeOperator {

	
	public final static String OPERATOR_NAME = ColumnOperator.class.getSimpleName();
	public final static OperationType OP = OperationType.COLUMN;
	
//	private int noOfHap;
	private int[] fixHaplotypeIndexArray;
	
	public ColumnOperator(HaplotypeModel haplotypeModel, int noOfHap, CoercionMode mode) {
		super(haplotypeModel, mode);

//		this.noOfHap = noOfHap;
//		this.noOfHap = haplotypeCount;
 
		fixHaplotypeIndexArray = new int[haplotypeCount];
		for (int i = 0; i < fixHaplotypeIndexArray.length; i++) {
			fixHaplotypeIndexArray[i] = i;
		}
	}

	@Override
	public String getPerformanceSuggestion() {

		return "";
	}

	@Override
	public double doOperation() throws OperatorFailedException {
		
		haplotypeModel.startAlignmentModelOperation();

		int siteIndex = getNextSiteIndex();
		
		for (int i = 0; i < haplotypeCount; i++) {
			Haplotype haplotype= haplotypeModel.getHaplotype(i);
			
			int oldState = haplotype.getState(siteIndex);
			char newChar = getNextDiffBase(oldState);
			haplotype.setCharAt(siteIndex, newChar);
		}

		haplotypeModel.setOperationRecord(OP, fixHaplotypeIndexArray, siteIndex);

		haplotypeModel.endAlignmentModelOperation();
		
		return 0.0;
		
		
		
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

	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

	@Override
	public OperationType getOperationType() {
		return OP;
	}

}
