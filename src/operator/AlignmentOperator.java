package operator;

import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;

public class AlignmentOperator extends AbstractCoercableOperator {

	public AlignmentOperator(CoercionMode mode) {
		super(mode);
		// TODO Auto-generated constructor stub
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
	public String getPerformanceSuggestion() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getOperatorName() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double doOperation() throws OperatorFailedException {
		// TODO Auto-generated method stub
		return 0;
	}


}
