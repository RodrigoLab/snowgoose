package srp.operator.haplotypes;

import java.util.Arrays;

import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;


public class HaplotypeSwapSectionOperator extends HaplotypeRecombinationOperator {

	//TODO: DEBUG: debug DEBUG?? to syntax tag here?
//	Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: -1
//    at srp.evolution.haplotypes.Haplotype.getChar(Haplotype.java:122)
//    at srp.operator.haplotypes.HaplotypeSwapSectionOperator.doOperation(HaplotypeSwapSectionOperator.java:59)
//    at dr.inference.operators.SimpleMCMCOperator.operate(SimpleMCMCOperator.java:172)
//    at dr.inference.markovchain.MarkovChain.runChain(MarkovChain.java:213)
//    at dr.inference.mcmc.MCMC.chain(MCMC.java:239)
//    at dr.inference.mcmc.MCMC.run(MCMC.java:195)
//    at srp.core.MainMCMCHaplotype.main(MainMCMCHaplotype.java:301)


	public final static String OPERATOR_NAME = HaplotypeSwapSectionOperator.class.getSimpleName();
//	public final static OperationType OP = OperationType.RECOMBINATION;

	
	
//	public AlignmentSwapBaseOperator(Parameter parameter, HaplotypeModel haplotypeModel, int index, CoercionMode mode) {
////		super(mode);
//	}


	public HaplotypeSwapSectionOperator(HaplotypeModel haplotypeModel) {
		super(haplotypeModel, (int) (haplotypeModel.getHaplotypeLength() * 0.01),
				CoercionMode.COERCION_OFF);
	}


	
	public HaplotypeSwapSectionOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
		super(haplotypeModel, length, mode);
	}

	@Override
	public double doOperation() throws OperatorFailedException {
		
		haplotypeModel.startAlignmentModelOperation();
		
		int[] twoHaplotypeIndex = new int[2];
		int[] twoPositionIndex = new int[2];
		
//		twoHaplotypeIndex[0] = getNextHapIndex();
//		twoHaplotypeIndex[1] = twoHaplotypeIndex[0];
		
		do{
			twoHaplotypeIndex[0] = getNextHapIndex();
			twoHaplotypeIndex[1] = getNextHapIndex();
		} while(twoHaplotypeIndex[0]==twoHaplotypeIndex[1]);
		
		do{
			twoPositionIndex[0] = getNextSiteIndex();
			twoPositionIndex[1]= twoPositionIndex[0]+ basesCount;
		}
		while(twoPositionIndex[1]>haplotypeLength);

		Haplotype h1 = haplotypeModel.getHaplotype(twoHaplotypeIndex[0]);
		Haplotype h2 = haplotypeModel.getHaplotype(twoHaplotypeIndex[1]);
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
