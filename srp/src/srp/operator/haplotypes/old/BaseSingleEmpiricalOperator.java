package srp.operator.haplotypes.old;

import srp.evolution.haplotypes.old.OldHaplotypeModel;
import dr.inference.operators.OperatorFailedException;

public class BaseSingleEmpiricalOperator extends AbstractBaseSingleOperator {

	public final static String OPERATOR_NAME = BaseSingleEmpiricalOperator.class.getSimpleName();

	public BaseSingleEmpiricalOperator(OldHaplotypeModel haplotypeModel,
			int nothing) {
		super(haplotypeModel);
	}

	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startHaplotypeOperation();

		int[] posChar = alignmentMapping.getNextBaseEmpirical();
		haplotypeModel.swapHaplotypeSingleBase(OP, posChar);

		haplotypeModel.endHaplotypeOperation();

		return 0.0;
	}

	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

}
