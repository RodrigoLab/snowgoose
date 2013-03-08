package srp.haplotypes.operator;

import srp.haplotypes.HaplotypeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;

public class AlignmentSwapBaseOperator extends AbstractCoercableOperator {

	int index;
	private HaplotypeModel haplotypeModel;
	
	public AlignmentSwapBaseOperator(Parameter parameter, HaplotypeModel haplotypeModel, int index, CoercionMode mode) {
		super(mode);
		// TODO Auto-generated constructor stub
	}

	
	public AlignmentSwapBaseOperator(HaplotypeModel haplotypeModel, int index, CoercionMode mode) {
		super(mode);
		this.index = index;
		this.haplotypeModel= haplotypeModel; 
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
		return "swapBaseOperator";
	}

	@Override
	public double doOperation() throws OperatorFailedException {
		haplotypeModel.swapBase();
		return 0;
	}

/*	
package operator;

import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;

import com.google.common.base.Strings;

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
	private ArrayList<String> alignment;

	public AlignmentOperator(Parameter variable, double scale, CoercionMode coercionOn, double d) {
		super(coercionOn);
		this.variable = variable;
		this.scale = scale;
		
	}
	
	public AlignmentOperator(Parameter variable, double scale) {

        this(variable, scale, CoercionMode.COERCION_ON, 1.0);
    }

	public AlignmentOperator(Parameter shortRead, ArrayList<String> haplotypes, int d) {
		this(shortRead, d);
		this.alignment = haplotypes;
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
		
		
		String t0 = alignment.get(0);
		StringBuilder sb = new StringBuilder(t0);
		
		String[] ACTG = new String[]{"A","C","G","T"};
		
		int l = (int) (Math.random()*100+500);
		t0 = StringUtils.repeat( ACTG[(int) (Math.random()*4)], l);
		
		String t1 = alignment.get(1);
		sb = new StringBuilder(t1);
		sb.setCharAt((int) d, 'P');
		l = (int) (Math.random()*100+500);
		t1 = StringUtils.repeat( ACTG[(int) (Math.random()*4)], l);
		
		System.out.println(d +"\t"+ t0.charAt(5) +"\t"+ t0.length() +"\t"+ t0);
		System.out.println(d +"\t"+ t1.charAt(5) +"\t"+ t1.length() +"\t"+ t1);
		
		alignment.set(0, t0);
		alignment.set(1, t1);
		
		
		
		
		return 0;
		
	}

	
	

}

	*/
}
