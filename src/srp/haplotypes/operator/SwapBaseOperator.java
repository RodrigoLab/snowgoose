package srp.haplotypes.operator;

import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.Operation;
import dr.evolution.alignment.Alignment;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

public class SwapBaseOperator extends SimpleMCMCOperator {

	public final static String OPERATOR_NAME = "SwapBaseOperator";
	public final static Operation OP = Operation.SWAPBASE;

	@Deprecated
	private int index;
	private HaplotypeModel haplotypeModel;
	
//	public AlignmentSwapBaseOperator(Parameter parameter, HaplotypeModel haplotypeModel, int index, CoercionMode mode) {
////		super(mode);
//	}

	
	public SwapBaseOperator(HaplotypeModel haplotypeModel, int nothing) {
//		super(mode);
		this.index = nothing;
		this.haplotypeModel= haplotypeModel; 
	}

	@Override
	public String getPerformanceSuggestion() {

//		System.err.println("getPero");
		return "getPerformanceSuggestion";
	}

	@Override
	public String getOperatorName() {

		return OPERATOR_NAME;
	}

	@Override
	public double doOperation() throws OperatorFailedException {
		
		haplotypeModel.startHaplotypeOperation();
			
		int hapIndex = MathUtils.nextInt( haplotypeModel.getHaplotypeCount());

		int[] posChar = haplotypeModel.getNextBase();
		int[] swapInfoArray = haplotypeModel.swapHaplotypeBase(hapIndex, posChar);
		
		haplotypeModel.storeOperationRecord(OP, swapInfoArray);
			
		haplotypeModel.endHaplotypeOperation();
		

		return 0.0;
	}

}
