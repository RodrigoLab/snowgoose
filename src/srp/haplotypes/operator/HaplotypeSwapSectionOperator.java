package srp.haplotypes.operator;

import org.apache.commons.math3.stat.StatUtils;

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


public class HaplotypeSwapSectionOperator extends AbstractSwapBasesOperator {

	
	public final static String OPERATOR_NAME = "SwapSectionOperator";
	public final static Operation OP = Operation.SWAPSECTION;

	
	
//	public AlignmentSwapBaseOperator(Parameter parameter, HaplotypeModel haplotypeModel, int index, CoercionMode mode) {
////		super(mode);
//	}

	
	public HaplotypeSwapSectionOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
		super(haplotypeModel, length, mode);

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
	public double getRawParameter() {
		return swapLength;
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
		int start = MathUtils.nextInt(haplotypeLength-swapLength+1);
		int end = start + swapLength; 


		String temp1 = oldS1.substring(start, end);
		String temp2 = oldS2.substring(start, end);

		String newS1 = oldS1.substring(0, start) + temp2 + oldS1.substring(end);
		String newS2 = oldS2.substring(0, start) + temp1 + oldS2.substring(end);
		
		h1.storeState();
		h1.setSequenceString(newS1);
		
		h2.storeState();
		h2.setSequenceString(newS2);
		
		int[] swapHaplotype = {hapIndex1, hapIndex2, start, end};

		haplotypeModel.storeOperationRecord(OP, swapHaplotype);
		haplotypeModel.endHaplotypeOperation();
		
		return 0.0;
		
		
		
	}

}
