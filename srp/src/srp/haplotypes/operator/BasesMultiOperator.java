package srp.haplotypes.operator;

import srp.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

public class BasesMultiOperator extends AbstractBasesMultiOperator {


	public final static String OPERATOR_NAME = BasesMultiOperator.class.getSimpleName();
//	public final static Operation OP = Operation.SWAPMULTI;

	

	
	public BasesMultiOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
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
			int[] posChar = alignmentMapping.getNextBase();
			allPosChars[0][posChar[0]] = posChar[1];
		}
		haplotypeModel.swapHaplotypeMultiBases(OP, allPosChars);

		haplotypeModel.endHaplotypeOperation();

		return 0.0;
	}

}
