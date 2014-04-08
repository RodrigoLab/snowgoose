package srp.operator.haplotypes.old;

import srp.evolution.haplotypes.old.OldHaplotypeModel;
import dr.inference.operators.OperatorFailedException;

public class BaseSingleOperator extends AbstractBaseSingleOperator {

	public final static String OPERATOR_NAME = BaseSingleOperator.class.getSimpleName();

	public BaseSingleOperator(OldHaplotypeModel haplotypeModel, int nothing) {
		super(haplotypeModel);
	}

	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startHaplotypeOperation();

		int[] posChar = alignmentMapping.getNextBase();
		haplotypeModel.swapHaplotypeSingleBase(OP, posChar);

		haplotypeModel.endHaplotypeOperation();

		return 0.0;
	}

	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

}
