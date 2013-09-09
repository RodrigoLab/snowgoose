package srp.haplotypes.operator;

import srp.haplotypes.HaplotypeModel;
import dr.inference.operators.OperatorFailedException;

public class SingleBaseEmpiricalOperator extends AbstractSingleBaseOperator {

	public final static String OPERATOR_NAME = SingleBaseEmpiricalOperator.class.getSimpleName();

	public SingleBaseEmpiricalOperator(HaplotypeModel haplotypeModel,
			int nothing) {
		super(haplotypeModel);
	}

	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startHaplotypeOperation();

		int[] posChar = alignmentMapping.getNextBaseEmpirical();
		haplotypeModel.swapHaplotypeSingleBase(OP, posChar);

		haplotypeModel.endHaplotypeOperation();

		return 0.0;
	}

}
