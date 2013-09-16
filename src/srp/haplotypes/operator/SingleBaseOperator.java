package srp.haplotypes.operator;

import srp.haplotypes.HaplotypeModel;
import dr.inference.operators.OperatorFailedException;

public class SingleBaseOperator extends AbstractSingleBaseOperator {

	public final static String OPERATOR_NAME = SingleBaseOperator.class.getSimpleName();

	public SingleBaseOperator(HaplotypeModel haplotypeModel, int nothing) {
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
