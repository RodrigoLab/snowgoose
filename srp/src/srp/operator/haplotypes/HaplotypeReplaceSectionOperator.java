package srp.operator.haplotypes;

import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;


public class HaplotypeReplaceSectionOperator extends HaplotypeRecombinationOperator {

	
	public final static String OPERATOR_NAME = HaplotypeReplaceSectionOperator.class.getSimpleName();
//	public final static OperationType OP = OperationType.RECOMBINATION;

	
	
//	public AlignmentSwapBaseOperator(Parameter parameter, HaplotypeModel haplotypeModel, int index, CoercionMode mode) {
////		super(mode);
//	}

	
	public HaplotypeReplaceSectionOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
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
//
		Haplotype h1 = haplotypeModel.getHaplotype(twoHaplotypeIndex[0]);
		Haplotype h2 = haplotypeModel.getHaplotype(twoHaplotypeIndex[1]);
		
		int maxStartingPos = haplotypeLength - basesCount;
		twoPositionIndex[0] = getNextSiteIndex(maxStartingPos);
		twoPositionIndex[1]= twoPositionIndex[0]+ basesCount;
		if(twoPositionIndex[1]>haplotypeLength){
			twoPositionIndex[1] = twoPositionIndex[0];
			twoPositionIndex[0] -= basesCount;
		}

		for (int i = twoPositionIndex[0]; i < twoPositionIndex[1]; i++) {
//			char c1 = h1.getChar(i);
			h1.setCharAt(i, h2.getChar(i));
//			h2.setCharAt(i, c1);
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
//FIXME: WRONG ratio? what is the right one??\
//		 System.out.println((1.0 / haplotypeCount));
//		 System.out.println((1.0 / Math.pow(4, basesCount))); 
		
		double xNewToOld = -Math.log(  Math.pow(4, basesCount) );
		double xOldToNew = -Math.log( (haplotypeCount-1) );
		double logq =  xNewToOld - xOldToNew  ;
		
		return logq;//TODO This might be wrong!! should not equal for the reverse step
		
		
		
	}



	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

}
