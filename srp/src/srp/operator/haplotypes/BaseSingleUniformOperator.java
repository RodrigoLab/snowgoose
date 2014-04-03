package srp.operator.haplotypes;

import srp.haplotypes.HaplotypeModel;
import dr.inference.operators.OperatorFailedException;

public class BaseSingleUniformOperator extends AbstractBaseSingleOperator {

	public final static String OPERATOR_NAME = BaseSingleUniformOperator.class.getSimpleName();

	public BaseSingleUniformOperator(HaplotypeModel haplotypeModel, int nothing) {
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
