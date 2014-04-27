package srp.operator.haplotypes;

import srp.evolution.OperationType;
import srp.haplotypes.Haplotype;
import srp.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

public class HaplotypeRecombinationOperator extends AbstractMultiOperator {

	public final static String OPERATOR_NAME = HaplotypeRecombinationOperator.class.getSimpleName();
	public final static OperationType OP = OperationType.RECOMBINATION;

	
	public HaplotypeRecombinationOperator(HaplotypeModel haplotypeModel, int nothing) {
		super(haplotypeModel, 0, CoercionMode.COERCION_OFF);
//		this.index = nothing;
//		this.haplotypeModel= haplotypeModel; 
//		this.haplotypeLength = this.haplotypeModel.getHaplotypeLength();
		
	}

	public HaplotypeRecombinationOperator(HaplotypeModel haplotypeModel,
			int length, CoercionMode mode) {
		super(haplotypeModel, length, mode);
	}

	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startAlignmentModelOperation();
		
		
		int[] twoHaplotypeIndex = new int[2];
		int[] twoPositionIndex = new int[2];
		
		
		twoHaplotypeIndex[0] = getNextHapIndex();
		twoHaplotypeIndex[1] = twoHaplotypeIndex[0];
		
		do{
			twoHaplotypeIndex[1] = getNextHapIndex();
		} while(twoHaplotypeIndex[0]==twoHaplotypeIndex[1]);

		Haplotype h1 = haplotypeModel.getHaplotype(twoHaplotypeIndex[0]);
		Haplotype h2 = haplotypeModel.getHaplotype(twoHaplotypeIndex[1]);
		
		String oldS1 = h1.getSequenceString();
		String oldS2 = h2.getSequenceString();
		
		do{
			twoPositionIndex[0] = getNextSiteIndex();
		}while(twoPositionIndex[0] ==0);

		twoPositionIndex[1]= haplotypeLength;


		String newS1 = oldS1.substring(0, twoPositionIndex[0]) + oldS2.substring(twoPositionIndex[0]);
		String newS2 = oldS2.substring(0, twoPositionIndex[0]) + oldS1.substring(twoPositionIndex[0]);
		
//		h1.storeState();
		h1.setSequenceString(newS1);
		
//		h2.storeState();
		h2.setSequenceString(newS2);
//		int[] swapHaplotype = {hapIndex1, hapIndex2, start, haplotypeLength};

		haplotypeModel.setOperationRecord(OP, twoHaplotypeIndex, twoPositionIndex);
		haplotypeModel.endAlignmentModelOperation();

		return 0.0;
	}


	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

	@Override
	public OperationType getOperationType() {
		return OP;
	}

}
