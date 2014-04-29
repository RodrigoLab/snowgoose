package srp.tdmcmc.evolution;

import srp.tdmcmc.SuperBeastTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.inference.model.Parameter;
import dr.math.MathUtils;


public class TransdimensionalTreeModel extends SuperBeastTreeModel{

	//
	// Public stuff
	//

	public TransdimensionalTreeModel(String name) {
		super(name);
	}
	public TransdimensionalTreeModel(Tree tree) {
		super(tree);
	}
	/**
	 * 
	 */
	private static final long serialVersionUID = -7300035033711283685L;

	public static final String TREE_MODEL = "treeModel";

	private static final boolean TEST_NODE_BOUNDS = false;

	int indexToPseudoRoot;//??

	private int trueInternalNodeCount;

	private int trueExternalNodeCount;

	private NodeRef trueRoot;

	private int trueNodeCount;
	public void InsertNode(Node n){
		
	}
	public void RemoveNode(Node n){
		
	}
	public void attachNode(Node n){
		
	}
	public void deattachNode(Node n){
		
	}

	public TransdimensionalTreeModel insertNode() {

//		System.out.println(this.toString());
		double logq = 0.0;

		NodeRef iGrandfather, iBrother;
		double heightFather;
		final int tipCount = this.getExternalNodeCount();
		final int internalNodeCount = this.getInternalNodeCount();

		final int nNodes = this.getNodeCount();
		NodeRef root = this.getRoot();

		int MAX_TRIES = 1000;
		int noSpareNode = externalNodeCount - trueExternalNodeCount;
		if(noSpareNode==0){
			noSpareNode = doubleTheTree();
		}
		trueNodeCount = externalNodeCount;
		Node parentNode;
		Node childNode;
		do {
			childNode = nodes[MathUtils.nextInt(trueNodeCount)];
		} while (root == childNode);
		System.err.println(nodeCount +"\t"+ nodes.length);
		
		parentNode = childNode.parent;

		
//		System.err.println(childNode.getNumber() + "\t" + childNode.getHeight()
//				+ "\t" + childNode.getChildCount());
//		System.err.println(parentNode.getNumber() + "\t"+parentNode.getChildCount());
//		System.err.println(parentNode==root);
//		System.err.println( "\t" + parentNode.getChild(0).getNumber() + "\t"
//				+ parentNode.getChild(1).getNumber());
		// FlexibleNode fn = (FlexibleNode) tree.getParent(childNode);
		// fn.getHeight();
		double u = MathUtils.nextDouble();
		double childHeight = childNode.getHeight();
		double parentHeight = parentNode.getHeight();
		double newExternalHeight = childHeight + u * (parentHeight - childHeight);
		double newInternalHeight = parentHeight - newExternalHeight;
		
//		newNode.setupHeightBounds();
//		newNode.setHeight(newNodeHeight);
//
//		System.err.println(parentNode.getNumber() + "\t"
//				+ parentNode.getHeight() + "\t" + parentNode.getChildCount());
//		System.err.println(newNode.getNumber() + "\t"
//				+ newNode.getHeight() + "\t" + newNode.getChildCount());
		

//		Node newExternalNode = createEmptyNode(trueExternalNodeCount);
		Node newExternalNode = new Node();
		newExternalNode.setNumber(trueExternalNodeCount);
		newExternalNode.heightParameter = new Parameter.Default(newExternalHeight);
		newExternalNode.taxon = new Taxon("XXXhap_"+trueExternalNodeCount);
		addVariable(newExternalNode.heightParameter);
//		newExternalNode.setHeight(newExternalHeight);
		nodes[trueExternalNodeCount] = newExternalNode;
		
//		Node newInternalNode = createEmptyNode(trueNodeCount+1);
		Node newInternalNode = new Node();
		newInternalNode.setNumber(trueNodeCount+1);
		newInternalNode.heightParameter = new Parameter.Default(newInternalHeight);
		newInternalNode.taxon = new Taxon("YYhap_"+trueNodeCount+1);
		addVariable(newInternalNode.heightParameter);
		nodes[trueNodeCount+1] = newInternalNode;
		
		System.out.println(nodes[trueExternalNodeCount].toString());
		System.out.println(nodes[trueNodeCount+1].toString());
//		newInternalNode.setHeight(newInternalHeight);
		
		trueExternalNodeCount++;
		trueInternalNodeCount++;
		trueNodeCount = trueExternalNodeCount+trueInternalNodeCount;
//		
//		Node newTip = new Node(this, parentNode);
//		newTip.removeChild(0);
//		newTip.removeChild(1);
//		newTip.number=0;
//		newNode.setHeight(newTipHeight);

//		for (int p = 0; p < parentNode.getChildCount(); p++) {
//			if (parentNode.getChild(p) == childNode) {
//				parentNode.removeChild(childNode);
//				System.out.println(p + "\t" + parentNode.getChildCount());
//			}
//		}
//		Node newNode = new Node(this, parentNode);
		
		
//		newNode.removeChild(0);
//		newNode.removeChild(1);
		newInternalNode.addChild(childNode);
		newInternalNode.addChild(newExternalNode);
		
		newExternalNode.parent = newInternalNode;
//		newNode.number=12;
//		newNode.parent = parentNode;
//		newTip.parent = newNode;
		
		parentNode.removeChild(childNode);
		parentNode.addChild(newInternalNode);

//		//
		System.err.println(parentNode.getNumber() + "\t"
				+ parentNode.getHeight() + "\t" + parentNode.getChildCount());

		System.err.println("child");
		System.err.println(childNode.getNumber() + "\t" + childNode.getHeight()
				+ "\t" + childNode.getChildCount());
		System.err.println("parent");
		System.err.println(parentNode.getNumber() + "\t"
				+ parentNode.getHeight() + "\t" + parentNode.getChildCount());
		System.err.println("\t" + parentNode.getChild(0).getNumber() + "\t"
				+ parentNode.getChild(1).getNumber());

		System.err.println("newNode");
		System.err.println(newInternalNode.getNumber() + "\t" + newInternalNode.getHeight()
				+ "\t" + newInternalNode.getChildCount());
		System.err.println("newTip");
		System.err.println(newExternalNode.getNumber() + "\t" + newExternalNode.getHeight()
				+ "\t" + newExternalNode.getChildCount());

		SimpleNode newRoot = new SimpleNode(this, root);
		SimpleTree st = new SimpleTree(newRoot);
		TransdimensionalTreeModel a = new TransdimensionalTreeModel(st);
		
		System.err.println(getInternalNodeCount() + "\t"
				+ getExternalNodeCount());
		System.err.println(st.getInternalNodeCount() + "\t"
				+ st.getExternalNodeCount());
		System.err.println(a.getInternalNodeCount() + "\t"
				+ a.getExternalNodeCount());

//		  +"\t"+ parentNode.getChild(0).getNumber() +"\t"+
//		 parentNode.getChild(1).getNumber());
//		// //
		// FlexibleNode fn = (FlexibleNode) treeModel.getRoot();
		// FlexibleTree ft = new FlexibleTree(fn);
		//
//		 System.err.println(st.toString());
//System.err.println("==========+");
//		root = st.getRoot();
//		 NodeRef node = root;
//			do {
//				node = Tree.Utils.postorderSuccessor(st, node);
////				System.err.println(node.getNumber());
//			} while (node != root);
		return a;
	}
	
    
	private int doubleTheTree() {
		throw new RuntimeException("Not yet implemented");
//		return 0;
	}

	public NodeRef getTrueRoot() {
        return trueRoot;
    }
	
	public int getTrueExternalNodeCount() {
		return trueExternalNodeCount;
	}
	public int getTrueInternalNodeCount() {
		return trueInternalNodeCount;
	}
	
	
	


}
