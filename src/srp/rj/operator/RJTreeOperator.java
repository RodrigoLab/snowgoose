package srp.rj.operator;

import java.util.Arrays;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.Operation;
import srp.haplotypes.SwapInfo;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeModel.Node;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

public class RJTreeOperator  extends AbstractTreeOperator {

	public final static Operation OP = Operation.TREE;
	
	protected HaplotypeModel haplotypeModel;

	protected AlignmentMapping alignmentMapping;

	private TreeModel tree;


	public final static String OPERATOR_NAME = RJTreeOperator.class.getSimpleName();

	public RJTreeOperator(HaplotypeModel haplotypeModel, TreeModel treeModel )  {

		this.haplotypeModel = haplotypeModel;
		this.tree = treeModel;
		
		
	}

//    public double getCoercableParameter() {
//        return Math.log(delta);
//    }
//
//    public void setCoercableParameter(double value) {
//        delta = Math.exp(value);
//    }


	/*
	*
	* @return the log-transformed hastings ratio
	*/
	@Override
	public double doOperation() throws OperatorFailedException {

		System.err.println(tree.toString());
		double logq = 0.0;

		NodeRef iGrandfather, iBrother;
		double heightFather;
		final int tipCount = tree.getExternalNodeCount();
		final int internalNodeCount = tree.getInternalNodeCount();

		final int nNodes = tree.getNodeCount();
		final NodeRef root = tree.getRoot();
		

		NodeRef i;

		int MAX_TRIES = 1000;

		
//		NodeRef parentNode = tree.getNode(MathUtils.nextInt(internalNodeCount));
	
//        for (int tries = 0; tries < MAX_TRIES; ++tries) {
			// get a random node whose father is not the root - otherwise
			// the operation is not possible
        	SimpleNode parentNode;
        	SimpleNode childNode;
			do {
				i=  tree.getNode(MathUtils.nextInt(nNodes));
			} while (root == i);
			childNode = new SimpleNode(tree, i);
			parentNode = new SimpleNode(tree, tree.getParent(i));
			System.err.println(childNode.getNumber()+"\t"+ childNode.getHeight() +"\t"+ 
			childNode.getChildCount());
			System.err.println(parentNode.getNumber()+"\t"+ parentNode.getHeight() +"\t"+ 
			parentNode.getChildCount() +"\t"+ parentNode.getChild(0).getNumber() +"\t"+ parentNode.getChild(1).getNumber());
//			FlexibleNode fn =  (FlexibleNode) tree.getParent(childNode);
//			fn.getHeight();
			double u = MathUtils.nextDouble();
			double childHeight = childNode.getHeight();
			double parentHeight = parentNode.getHeight();
			double newTipHeight = childHeight + u*(parentHeight-childHeight);
			double newNodeHeight = parentHeight - newTipHeight;
			
			SimpleNode newNode = new SimpleNode();
			newNode.setHeight(newNodeHeight);
			
			SimpleNode newTip = new SimpleNode ();
			newNode.setHeight(newTipHeight);
			
			for (int p = 0; p < parentNode.getChildCount(); p++) {
				if(parentNode.getChild(p)==childNode){
					parentNode.removeChild(childNode);
					System.out.println(p +"\t"+ parentNode.getChildCount());
				}
			}
//			parentNode.removeChild(childNode.getNumber());
			parentNode.addChild(newNode);
			
			newNode.addChild(childNode);
			newNode.addChild(newTip);
			
			System.err.println(parentNode.getNumber()+"\t"+ parentNode.getHeight() +"\t"+ 
					parentNode.getChildCount() ); 
//					+"\t"+ parentNode.getChild(0).getNumber() +"\t"+ parentNode.getChild(1).getNumber());
//					
			FlexibleNode fn = (FlexibleNode) tree.getRoot();
			FlexibleTree ft = new FlexibleTree(fn);
			
			System.err.println(ft.toString());
           // int childIndex = (MathUtils.nextDouble() >= 0.5 ? 1 : 0);
           // int otherChildIndex = 1 - childIndex;
           // NodeRef iOtherChild = tree.getChild(i, otherChildIndex);
//
//           NodeRef iFather = tree.getParent(i);
//           iGrandfather = tree.getParent(iFather);
//           iBrother = getOtherChild(tree, iFather, i);
//           heightFather = tree.getNodeHeight(iFather);
//
//           // NodeRef newChild = getRandomNode(possibleChilds, iFather);
//           NodeRef newChild = tree.getNode(MathUtils.nextInt(nNodes));

//           if (tree.getNodeHeight(newChild) < heightFather
//                 && root != newChild
//                 && tree.getNodeHeight(tree.getParent(newChild)) > heightFather
//                 && newChild != iFather
//                 && tree.getParent(newChild) != iFather) {
//              NodeRef newGrandfather = tree.getParent(newChild);

//              tree.beginTreeEdit();
//
//              // prune
//              tree.removeChild(iFather, iBrother);
//              tree.removeChild(iGrandfather, iFather);
//              tree.addChild(iGrandfather, iBrother);
//
//              // reattach
//              tree.removeChild(newGrandfather, newChild);
//              tree.addChild(iFather, newChild);
//              tree.addChild(newGrandfather, iFather);
//
//              // ****************************************************
//
//              tree.endTreeEdit();
//
//              tree.pushTreeChangedEvent(i);
//
//              assert tree.getExternalNodeCount() == tipCount;

//        }
		System.exit(-1);
		return logq;
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
