package srp.haplotypes.operator;

import java.util.Arrays;
import java.util.HashMap;

import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.Operation;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class SwapBasesUniformOperator extends AbstractHaplotypeOperator{


	public final static String OPERATOR_NAME = SwapBasesUniformOperator.class.getSimpleName();
	public final static Operation OP = Operation.SWAPMULTI;

	private int[][] allPosChars; 
	private int[] allNewChars;
	
	public SwapBasesUniformOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
		super(haplotypeModel, length, mode);
		allPosChars = new int[2][haplotypeLength];
		allNewChars = new int[haplotypeLength];
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

		int hapIndex = MathUtils.nextInt(haplotypeCount);
		
		for (int i = 0; i < allPosChars.length; i++) {
			Arrays.fill(allPosChars[i], -1);
		}

		for (int i = 0; i < swapLength; i++) {
			int[] posChar = haplotypeModel.getNextBaseUniform();
			allPosChars[0][posChar[0]] = posChar[1];
		}
		
		haplotypeModel.swapHaplotypeMultiBases(hapIndex, allPosChars[0], allPosChars[1]);

		haplotypeModel.storeOperationRecord(OP, hapIndex, allPosChars);
		
		haplotypeModel.endHaplotypeOperation();

		return 0.0;
	}
	
}