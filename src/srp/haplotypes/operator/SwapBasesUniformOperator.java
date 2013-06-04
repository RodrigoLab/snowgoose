package srp.haplotypes.operator;

import java.util.Arrays;
import java.util.HashMap;

import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.Operation;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class SwapBasesUniformOperator extends AbstractSwapBasesOperator{


	public final static String OPERATOR_NAME = SwapBasesUniformOperator.class.getSimpleName();
	public final static Operation OP = Operation.SWAPMULTI;

	
	public SwapBasesUniformOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
		super(haplotypeModel, length, mode);

	}

    @Override
	public double getRawParameter() {
		// 
//		System.err.println("getRaw");
		return swapLength;
	}

	@Override
	public String getPerformanceSuggestion() {

		return "getPerformanceSuggestion";
	}
	@Override
	public String getOperatorName() {

		return OPERATOR_NAME;
	}

	
	
	@Override
	public double doOperation() throws OperatorFailedException {
		
		haplotypeModel.startHaplotypeOperation();
		
		for (int i = 0; i < allPosChars.length; i++) {
			Arrays.fill(allPosChars[i], -1);
		}
		for (int i = 0; i < swapLength; i++) {
			int[] posChar = alignmentMapping.getNextBaseUniform();
			allPosChars[0][posChar[0]] = posChar[1];
		}
		
		int hapIndex = MathUtils.nextInt(haplotypeCount);
		haplotypeModel.swapHaplotypeMultiBases(hapIndex, allPosChars[0], allPosChars[1]);

		haplotypeModel.storeOperationRecord(OP, hapIndex, allPosChars);
		
		haplotypeModel.endHaplotypeOperation();

		return 0.0;
	}
	
}