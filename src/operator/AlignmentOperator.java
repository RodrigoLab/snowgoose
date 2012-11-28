package operator;

import dr.evolution.alignment.SimpleAlignment;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.OperatorUtils;
import dr.inference.operators.SimpleMCMCOperator;

public class AlignmentOperator extends AbstractCoercableOperator {

	private double scale;
	private Parameter variable;
	private SimpleAlignment alignment;

	public AlignmentOperator(Parameter variable, double scale, CoercionMode coercionOn, double d) {
		super(coercionOn);
		this.variable = variable;
		this.scale = scale;
		
	}
	
	public AlignmentOperator(Parameter variable, double scale) {

        this(variable, scale, CoercionMode.COERCION_ON, 1.0);
    }

	public AlignmentOperator(Parameter shortRead, SimpleAlignment alignment, int d) {
		this(shortRead, d);
		this.alignment = alignment;
	}


	@Override
	public double getCoercableParameter() {
//        return Math.log(1.0 / scaleFactor - 1.0);
		return scale;

	}

	@Override
	public void setCoercableParameter(double value) {
//		scaleFactor = 1.0 / (Math.exp(value) + 1.0);
		scale = value;
		
	}

	@Override
	public double getRawParameter() {
		return scale;
	}

	@Override
	public String getPerformanceSuggestion() {
//        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
//        double targetProb = getTargetAcceptanceProbability();
//        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
//        double sf = OperatorUtils.optimizeScaleFactor(scaleFactor, prob, targetProb);
//        if (prob < getMinimumGoodAcceptanceLevel()) {
//            return "Try setting scaleFactor to about " + formatter.format(sf);
//        } else if (prob > getMaximumGoodAcceptanceLevel()) {
//            return "Try setting scaleFactor to about " + formatter.format(sf);
//        } else return "";
		return "NOT YET IMPLEMETED";
	}

	@Override
	public String getOperatorName() {
//		return "scale(" + variable.getVariableName() + ")";
		
		return "AlignmentOperator....";
	}

	@Override
	public double doOperation() throws OperatorFailedException {
		double d = variable.getValue(0);
		
		variable.setValue(0, d+1);
		System.out.println(d +"\t"+ variable.getValue(0));
		return 0;
	}


}
