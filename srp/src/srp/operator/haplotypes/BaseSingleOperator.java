package srp.operator.haplotypes;

import srp.evolution.OperationType;
import srp.haplotypes.Haplotype;
import srp.haplotypes.HaplotypeModel;
import dr.inference.operators.OperatorFailedException;

public class BaseSingleOperator extends AbstractSingleOperator {

	public static final String OPERATOR_NAME = BaseSingleOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.SINGLE;
	
	public BaseSingleOperator(HaplotypeModel haplotypeModel) {
		super(haplotypeModel);
	}

	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startAlignmentModelOperation();

		int hapIndex = getNextHapIndex();
		Haplotype haplotype = haplotypeModel.getHaplotype(hapIndex);
		
		int siteIndex = getNextSiteIndex();
		int oldState = haplotype.getState(siteIndex);
		char newChar = getNextDiffBase(oldState);
		
		haplotype.setCharAt(siteIndex, newChar);
		haplotypeModel.setOperationRecord(OP, hapIndex, siteIndex);

		haplotypeModel.endAlignmentModelOperation();

		return 0.0;
	}

	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}




}
