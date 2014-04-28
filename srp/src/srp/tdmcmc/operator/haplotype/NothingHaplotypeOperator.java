package srp.tdmcmc.operator.haplotype;

import java.util.Arrays;

import srp.evolution.haplotypes.old.OldHapOperation;
import srp.evolution.haplotypes.old.OldHapSwapInfo;
import srp.evolution.haplotypes.old.OldHaplotypeModel;
import srp.evolution.shortreads.AlignmentMapping;
import srp.tdmcmc.evolution.TransdimensionalHaplotypeModel;
import srp.tdmcmc.evolution.TransdimensionalTreeModel;
import srp.tdmcmc.operator.AbstractTransdimensionalOperator;
import dr.evolution.datatype.Nucleotides;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class NothingHaplotypeOperator extends AbstractTransdimensionalOperator {

	public final static String OPERATOR_NAME = NothingHaplotypeOperator.class.getSimpleName();

	private TransdimensionalTreeModel tree;
	private TransdimensionalHaplotypeModel haplotypeModel;

	public NothingHaplotypeOperator(TransdimensionalHaplotypeModel haplotypeModel, double param, CoercionMode mode )  {
		super(mode);

	
	}

	@Override
	public double doOperation() throws OperatorFailedException {
		// TODO Auto-generated method stub
		return 0;
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
		return OPERATOR_NAME;
	}
	
}
