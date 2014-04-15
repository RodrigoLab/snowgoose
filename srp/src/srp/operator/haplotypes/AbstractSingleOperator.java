package srp.operator.haplotypes;

import srp.evolution.OperationType;
import srp.haplotypes.HaplotypeModel;

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
	
}