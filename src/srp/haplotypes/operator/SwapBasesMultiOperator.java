package srp.haplotypes.operator;

import java.util.Arrays;

import srp.haplotypes.Haplotype;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.Operation;
import dr.evolution.alignment.Alignment;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

public class SwapBasesMultiOperator extends AbstractSwapBasesOperator {


	public final static String OPERATOR_NAME = SwapBasesMultiOperator.class.getSimpleName();
//	public final static Operation OP = Operation.SWAPMULTI;

	

	
	public SwapBasesMultiOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
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
