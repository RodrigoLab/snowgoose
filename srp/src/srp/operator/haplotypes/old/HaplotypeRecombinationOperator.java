package srp.operator.haplotypes.old;

import srp.evolution.haplotypes.old.OldHapOperation;
import srp.evolution.haplotypes.old.OldHaplotype;
import srp.evolution.haplotypes.old.OldHaplotypeModel;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

public class HaplotypeRecombinationOperator extends SimpleMCMCOperator {

	public final static String OPERATOR_NAME = HaplotypeRecombinationOperator.class.getSimpleName();
	public final static OldHapOperation OP = OldHapOperation.SWAPSECTION;

	@Deprecated
	private int index;
	private int haplotypeLength;
	private OldHaplotypeModel haplotypeModel;
	
//	public AlignmentSwapBaseOperator(Parameter parameter, HaplotypeModel haplotypeModel, int index, CoercionMode mode) {
////		super(mode);
//	}

	
	public HaplotypeRecombinationOperator(OldHaplotypeModel haplotypeModel, int nothing) {
//		super(mode);
		this.index = nothing;
		this.haplotypeModel= haplotypeModel; 
		this.haplotypeLength = this.haplotypeModel.getHaplotypeLength();
		
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
			
		int hapIndex1 = MathUtils.nextInt( haplotypeModel.getHaplotypeCount());
		int hapIndex2 = hapIndex1;
		
		do{
			hapIndex2 = MathUtils.nextInt( haplotypeModel.getHaplotypeCount());
		} while(hapIndex1==hapIndex2);

		OldHaplotype h1 = haplotypeModel.getHaplotype(hapIndex1);
		OldHaplotype h2 = haplotypeModel.getHaplotype(hapIndex2);
		
		String oldS1 = h1.getSequenceString();
		String oldS2 = h2.getSequenceString();
		int start = MathUtils.nextInt(haplotypeLength);
		do{
			start = MathUtils.nextInt(haplotypeLength);
		}while(start==0);

		String newS1 = oldS1.substring(0, start) + oldS2.substring(start);
		String newS2 = oldS2.substring(0, start) + oldS1.substring(start);
		
		h1.storeState();
		h1.setSequenceString(newS1);
		
		h2.storeState();
		h2.setSequenceString(newS2);
		int[] swapHaplotype = {hapIndex1, hapIndex2, start, haplotypeLength};

		haplotypeModel.storeOperationRecord(OP, swapHaplotype);
		haplotypeModel.endHaplotypeOperation();

		return 0.0;
	}

}
