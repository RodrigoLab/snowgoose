package srp.operator.haplotypes.old;

import srp.evolution.haplotypes.old.OldHaplotypeModel;
import dr.inference.operators.OperatorFailedException;

public class BaseSingleUniformOperator extends AbstractBaseSingleOperator {

	public final static String OPERATOR_NAME = BaseSingleUniformOperator.class.getSimpleName();

	public BaseSingleUniformOperator(OldHaplotypeModel haplotypeModel, int nothing) {
		super(haplotypeModel);

	}

	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startHaplotypeOperation();

		int[] posChar = alignmentMapping.getNextBaseUniform();
		
		haplotypeModel.swapHaplotypeSingleBase(OP, posChar);

		haplotypeModel.endHaplotypeOperation();

		return 0.0;
	}

	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

}
