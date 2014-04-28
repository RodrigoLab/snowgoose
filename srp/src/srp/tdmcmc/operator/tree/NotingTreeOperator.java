package srp.tdmcmc.operator.tree;

import srp.evolution.haplotypes.old.OldHapOperation;
import srp.evolution.haplotypes.old.OldHaplotypeModel;
import srp.evolution.shortreads.AlignmentMapping;
import srp.tdmcmc.evolution.TransdimensionalHaplotypeModel;
import srp.tdmcmc.evolution.TransdimensionalTreeModel;
import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class NotingTreeOperator  extends AbstractTreeOperator {

	public final static String OPERATOR_NAME = MoveRootTreeOperator.class.getSimpleName();

	private TransdimensionalTreeModel tree;
	private TransdimensionalHaplotypeModel haplotypeModel;


	public NotingTreeOperator(TransdimensionalHaplotypeModel haplotypeModel, 
			TransdimensionalTreeModel treeModel )  {

		this.haplotypeModel = haplotypeModel;
		this.tree = treeModel;
		
		
	}


	/*
	*
	* @return the log-transformed hastings ratio
	*/
	@Override
	public double doOperation() throws OperatorFailedException {
		return 0;
	}
	
	

	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}
//
//    @Override
//	public double getMinimumAcceptanceLevel() {
//        return targetAcceptanceProb-0.2;
//    }
//
//    @Override
//	public double getMaximumAcceptanceLevel() {
//        return targetAcceptanceProb + 0.2;
//    }
//
//    @Override
//	public double getMinimumGoodAcceptanceLevel() {
//        return targetAcceptanceProb- 0.1;
//    }
//
//    @Override
//	public double getMaximumGoodAcceptanceLevel() {
//        return targetAcceptanceProb+ 0.1;
//    }
//
//    private double targetAcceptanceProb = 0.5;
//
//
//	

	@Override
	public String getPerformanceSuggestion() {

		return null;
	}
}
