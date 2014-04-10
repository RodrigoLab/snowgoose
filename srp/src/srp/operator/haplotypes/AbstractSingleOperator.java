package srp.operator.haplotypes;

import srp.evolution.OperationType;
import srp.haplotypes.HaplotypeModel;
import dr.math.MathUtils;

public abstract class AbstractSingleOperator extends AbstractHaplotypeOperator {

	public static final OperationType OP = OperationType.SINGLE;
	
	public AbstractSingleOperator(HaplotypeModel haplotypeModel) {
		super(haplotypeModel);
		
	}

	public char getNextDiffBase(int oldState){
		int i = oldState;
		do {
			i = MathUtils.nextInt(DIMENSION);
		} while (i == oldState);

		return DNA_CHARS[i];
	}

	public OperationType getOperationType(){
		return OP;
	}


	@Override
	public String getPerformanceSuggestion() {
		return "";
	}
	
}