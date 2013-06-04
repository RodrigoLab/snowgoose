package srp.haplotypes.operator;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.Operation;
import dr.evolution.alignment.Alignment;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;


public class SingleBaseOperator extends SimpleMCMCOperator {

	public final static String OPERATOR_NAME = SingleBaseOperator.class.getSimpleName();
	public final static Operation OP = Operation.SWAPSINGLE;

	@Deprecated
	private int index;
	private HaplotypeModel haplotypeModel;
	private int haplotypeCount;
	protected AlignmentMapping alignmentMapping;
//	public AlignmentSwapBaseOperator(Parameter parameter, HaplotypeModel haplotypeModel, int index, CoercionMode mode) {
////		super(mode);
//	}

	
	public SingleBaseOperator(HaplotypeModel haplotypeModel, int nothing) {
//		super(mode);
		this.index = nothing;
		this.haplotypeModel= haplotypeModel; 
		haplotypeCount = this.haplotypeModel.getHaplotypeCount();
		alignmentMapping = this.haplotypeModel.getAlignmentMapping();
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
			
		int hapIndex = MathUtils.nextInt( haplotypeCount);

		int[] posChar = alignmentMapping.getNextBase();
		int[] swapInfoArray = haplotypeModel.swapHaplotypeSingleBase(hapIndex, posChar);
		
		haplotypeModel.storeOperationRecord(OP, swapInfoArray);
			
		haplotypeModel.endHaplotypeOperation();
		

		return 0.0;
	}

}
