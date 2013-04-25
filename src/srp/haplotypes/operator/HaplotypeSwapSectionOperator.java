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


public class HaplotypeSwapSectionOperator extends AbstractCoercableOperator {

	
	public final static String OPERATOR_NAME = "SwapSectionOperator";
	public final static Operation OP = Operation.SWAPSECTION;

	
	private int swapFragmentLength;
	private int hapLength;
	private HaplotypeModel haplotypeModel;
	
	
//	public AlignmentSwapBaseOperator(Parameter parameter, HaplotypeModel haplotypeModel, int index, CoercionMode mode) {
////		super(mode);
//	}

	
	public HaplotypeSwapSectionOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
		super(mode);
		this.swapFragmentLength = length;
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
	public double getCoercableParameter() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setCoercableParameter(double value) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double getRawParameter() {
		// TODO Auto-generated method stub
		return 0;
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
		int pos = MathUtils.nextInt(hapLength-swapFragmentLength+1);
		int end = pos + swapFragmentLength; 


		String temp1 = oldS1.substring(pos, end);
		String temp2 = oldS2.substring(pos, end);

		String newS1 = oldS1.substring(0, pos) + temp2 + oldS1.substring(end);
		String newS2 = oldS2.substring(0, pos) + temp1 + oldS2.substring(end);
		
		h1.storeState();
		h1.setSequenceString(newS1);
		
		h2.storeState();
		h2.setSequenceString(newS2);
		
		int[] swapHaplotype = {hapIndex1, hapIndex2};

		haplotypeModel.storeOperationRecord(OP, swapHaplotype);
		haplotypeModel.endHaplotypeOperation();
		
		return 0.0;
		
		
		
	}

}
