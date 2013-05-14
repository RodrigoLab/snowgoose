package srp.haplotypes.operator;

import java.util.ArrayList;
import java.util.Arrays;

import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.Operation;
import dr.evolution.alignment.Alignment;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

public class SingleBaseUniformOperator extends SimpleMCMCOperator {

	public final static String OPERATOR_NAME = SingleBaseUniformOperator.class.getSimpleName();
	public final static Operation OP = Operation.SWAPSINGLE;

	@Deprecated
	private int index;
	private HaplotypeModel haplotypeModel;
	private int haplotypeCount;
	

	
	public SingleBaseUniformOperator(HaplotypeModel haplotypeModel, int nothing) {
//		super(mode);
		this.index = nothing;
		this.haplotypeModel= haplotypeModel; 
		haplotypeCount = this.haplotypeModel.getHaplotypeCount();
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
			
		int hapIndex = MathUtils.nextInt(haplotypeCount);

		int[] posChar = haplotypeModel.getNextBaseUniform();
		int[] swapInfoArray = haplotypeModel.swapHaplotypeSingleBase(hapIndex, posChar);
		
		haplotypeModel.storeOperationRecord(OP, swapInfoArray);
			
		haplotypeModel.endHaplotypeOperation();
		

		return 0.0;
	}


}
