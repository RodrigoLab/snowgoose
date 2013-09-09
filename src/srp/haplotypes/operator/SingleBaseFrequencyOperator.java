package srp.haplotypes.operator;

import srp.haplotypes.HaplotypeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.OperatorFailedException;

public class SingleBaseFrequencyOperator extends AbstractSingleBaseOperator {

	public final static String OPERATOR_NAME = SingleBaseFrequencyOperator.class.getSimpleName();
	private Parameter frequency;

	public SingleBaseFrequencyOperator(HaplotypeModel haplotypeModel,
			Parameter freqs) {

		super(haplotypeModel);
		this.frequency = freqs;
	}

	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startHaplotypeOperation();

		int[] posChar = alignmentMapping.getNextBaseFrequency(frequency);
		haplotypeModel.swapHaplotypeSingleBase(OP, posChar);

		haplotypeModel.endHaplotypeOperation();

		return 0.0;
	}

}
