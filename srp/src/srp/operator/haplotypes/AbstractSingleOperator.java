package srp.operator.haplotypes;

import srp.evolution.OperationType;
import srp.evolution.haplotypes.HaplotypeModel;

public abstract class AbstractSingleOperator extends AbstractHaplotypeOperator {

	public static final OperationType OP = OperationType.SINGLE;
	
	public AbstractSingleOperator(HaplotypeModel haplotypeModel) {
		super(haplotypeModel);
		
	}

	public OperationType getOperationType(){
		return OP;
	}


	@Override
	public String getPerformanceSuggestion() {
		return "";
	}

	@Override
	public double getCoercableParameter() {
		return 0;
	}

	@Override
	public void setCoercableParameter(double value) {
	}

	@Override
	public double getRawParameter() {
		return 0;
	}

	
}