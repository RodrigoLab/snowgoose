package srp.haplotypes.operator;

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

public class HaplotypeRecombinationOperator extends SimpleMCMCOperator {

	public final static String OPERATOR_NAME = "HaplotypeRecombinationOperator";
	public final static Operation OP = Operation.RECOMB;

	@Deprecated
	private int index;
	private int hapLength;
	private HaplotypeModel haplotypeModel;
	
//	public AlignmentSwapBaseOperator(Parameter parameter, HaplotypeModel haplotypeModel, int index, CoercionMode mode) {
////		super(mode);
//	}

	
	public HaplotypeRecombinationOperator(HaplotypeModel haplotypeModel, int nothing) {
//		super(mode);
		this.index = nothing;
		this.haplotypeModel= haplotypeModel; 
		this.hapLength = this.haplotypeModel.getHaplotypeLength();
		
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
			
		int hapIndex1 = MathUtils.nextInt( haplotypeModel.getHaplotypeCount());
		int hapIndex2 = hapIndex1;
		
		do{
			hapIndex2 = MathUtils.nextInt( haplotypeModel.getHaplotypeCount());
		} while(hapIndex1==hapIndex2);

		Haplotype h1 = haplotypeModel.getHaplotype(hapIndex1);
		Haplotype h2 = haplotypeModel.getHaplotype(hapIndex2);
		
		String oldS1 = h1.getSequenceString();
		String oldS2 = h2.getSequenceString();
		int start = MathUtils.nextInt(hapLength);
		do{
			start = MathUtils.nextInt(hapLength);
		}while(start==0);

		String newS1 = oldS1.substring(0, start) + oldS2.substring(start);
		String newS2 = oldS2.substring(0, start) + oldS1.substring(start);
		
		h1.storeState();
		h1.setSequenceString(newS1);
		
		h2.storeState();
		h2.setSequenceString(newS2);
		
		int[] swapHaplotype = {hapIndex1, hapIndex2, start, hapLength};

		haplotypeModel.storeOperationRecord(OP, swapHaplotype);
		haplotypeModel.endHaplotypeOperation();

		return 0.0;
	}

}
