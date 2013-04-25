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

public class SwapMultiBasesOperator extends AbstractCoercableOperator {


	public final static String OPERATOR_NAME = "SwapMultiBasesOperator";
	public final static Operation OP = Operation.SWAPMULTI;


	int swapNBases;
	private HaplotypeModel haplotypeModel;
	
	public SwapMultiBasesOperator(Parameter parameter, HaplotypeModel haplotypeModel, int swapNBase, CoercionMode mode) {
		super(mode);
		// TODO Auto-generated constructor stub
	}

	
	public SwapMultiBasesOperator(HaplotypeModel haplotypeModel, int nBases, CoercionMode mode) {
		super(mode);
		if(nBases<1){
			nBases=1;
		}
		this.swapNBases =  nBases;
		this.haplotypeModel= haplotypeModel; 
	}

	@Override
	public double getCoercableParameter() {

		System.out.println("getCoercableP: "+swapNBases);
		return Math.log(swapNBases);
	}

	@Override
	public void setCoercableParameter(double value) {
		swapNBases = (int) (Math.exp(value));
//		if (swapNBase == 0)
//			swapNBase = 1;
		System.out.println("setCoer\t"+ swapNBases +"\t"+ getAcceptanceProbability());
	}

	
	@Override
	public double getRawParameter() {
		// 
//		System.err.println("getRaw");
		return swapNBases;
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
//		haplotypeModel.swapMultiBase(swapNBase);

		haplotypeModel.startHaplotypeOperation();

		int hapIndex = MathUtils.nextInt( haplotypeModel.getHaplotypeCount());
		haplotypeModel.storeOperationRecord(OP, null);
		for (int i = 0; i < swapNBases; i++) {

			int[] posChar = haplotypeModel.getNextBase();
			int[] swapInfoArray = haplotypeModel.swapHaplotypeBase(hapIndex, posChar);

			haplotypeModel.storeOperationRecord(OP, swapInfoArray);
		}

		haplotypeModel.endHaplotypeOperation();
		
		

		
		
		return 0.0;
	}

}