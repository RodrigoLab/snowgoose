package srp.operator.haplotypes;

import srp.evolution.OperationType;
import srp.evolution.haplotypes.old.OldHaplotype;
import srp.evolution.haplotypes.old.OldHaplotypeModel;
import srp.haplotypes.Haplotype;
import srp.haplotypes.HaplotypeModel;
import srp.operator.haplotypes.old.AbstractBasesMultiOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;


public class HaplotypeSwapSectionOperator extends HaplotypeRecombinationOperator {

	
	public final static String OPERATOR_NAME = HaplotypeSwapSectionOperator.class.getSimpleName();
//	public final static OperationType OP = OperationType.RECOMBINATION;

	
	
//	public AlignmentSwapBaseOperator(Parameter parameter, HaplotypeModel haplotypeModel, int index, CoercionMode mode) {
////		super(mode);
//	}

	
	public HaplotypeSwapSectionOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
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
//FIXME fix index -negative index
		Haplotype h1 = haplotypeModel.getHaplotype(twoHaplotypeIndex[0]);
		Haplotype h2 = haplotypeModel.getHaplotype(twoHaplotypeIndex[1]);
		
		twoPositionIndex[0] = getNextSiteIndex();
		twoPositionIndex[1]= twoPositionIndex[0]+ basesCount;
		if(twoPositionIndex[1]>haplotypeLength){
			twoPositionIndex[1] = twoPositionIndex[0];
			twoPositionIndex[0] -= basesCount;
		}

		for (int i = twoPositionIndex[0]; i < twoPositionIndex[1]; i++) {
			char c1 = h1.getChar(i);
			h1.setCharAt(i, h2.getChar(i));
			h2.setCharAt(i, c1);
		}

//		int start = twoPositionIndex[0];
//		int end = twoPositionIndex[1];
//		String oldS1 = h1.getSequenceString();
//		String oldS2 = h2.getSequenceString();
//		String temp1 = oldS1.substring(start, end);
//		String temp2 = oldS2.substring(start, end);
//		String newS1 = oldS1.substring(0, start) + temp2 + oldS1.substring(end);
//		String newS2 = oldS2.substring(0, start) + temp1 + oldS2.substring(end);
//		h1.setSequenceString(newS1);
//		h2.setSequenceString(newS2);

		haplotypeModel.setOperationRecord(OP, twoHaplotypeIndex, twoPositionIndex);
		haplotypeModel.endAlignmentModelOperation();

		return 0.0;
		
		
		
	}



	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

}
