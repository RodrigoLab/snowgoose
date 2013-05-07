package srp.haplotypes.operator;

import java.util.ArrayList;

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

public class SwapBasesUniformOperator extends AbstractHaplotypeOperator{


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
//		haplotypeModel.swapMultiBase(swapNBase);

		haplotypeModel.startHaplotypeOperation();

		int hapIndex = MathUtils.nextInt( haplotypeModel.getHaplotypeCount());
		haplotypeModel.storeOperationRecord(OP, null);
		for (int i = 0; i < swapLength; i++) {

			int[] posChar = haplotypeModel.getNextBaseUniform();
			int[] swapInfoArray = haplotypeModel.swapHaplotypeBase(hapIndex, posChar);

			haplotypeModel.storeOperationRecord(OP, swapInfoArray);
		}

		haplotypeModel.endHaplotypeOperation();
		
		

		
		
		return 0.0;
	}
	
}