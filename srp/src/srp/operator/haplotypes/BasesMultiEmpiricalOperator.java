package srp.operator.haplotypes;

import srp.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

public class BasesMultiEmpiricalOperator extends AbstractBasesMultiOperator {


	public final static String OPERATOR_NAME = BasesMultiEmpiricalOperator.class.getSimpleName();
	
	public BasesMultiEmpiricalOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
		super(haplotypeModel, length, mode);
		
		
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
		resetAllPosChars();

		for (int i = 0; i < swapLength; i++) {
			int[] posChar = alignmentMapping.getNextBaseEmpirical();
			allPosChars[0][posChar[0]] = posChar[1];
		}
		haplotypeModel.swapHaplotypeMultiBases(OP, allPosChars);

		haplotypeModel.endHaplotypeOperation();

		return 0.0;
	}

}