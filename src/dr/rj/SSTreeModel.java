package dr.rj;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.swing.text.TabableView;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.google.common.util.concurrent.SimpleTimeLimiter;

import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.MutableTreeListener;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.inference.model.AbstractModel;
import dr.inference.model.Bounds;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.inference.model.Variable;
import dr.math.MathUtils;
import dr.math.distributions.ReflectedNormalDistribution;
import dr.util.Attributable;


public class SSTreeModel extends AbstractModel implements MutableTree {

	//
	// Public stuff
	//
	public static final String TREE_MODEL = "treeModel";

	private static final boolean TEST_NODE_BOUNDS = false;
	//
	// SSTreeModel
	//
	private static final int SCALE_TAXA_FACTOR = 2;
	private int initExternalNode = 16;
//	private int initInternalNode = initExternalNode-1;

	private NodeRef trueRoot;
	private int trueExternalNodeCount;
	private int trueInternalNodeCount;
	private int trueNodeCount;
	

	public SSTreeModel(String name) {
		super(name);
		nodeCount = 0;
		externalNodeCount = 0;
		internalNodeCount = 0;
	}

	public SSTreeModel(Tree tree) {
		this(TREE_MODEL, tree, false);
	}

	public SSTreeModel(String id, Tree tree) {

		this(TREE_MODEL, tree, false);
		setId(id);
	}

	/*
	 * New constructor that copies the attributes of Tree tree into the new
	 * TreeModel Useful for constructing a TreeModel from a NEXUS file entry
	 */

	public SSTreeModel(String name, Tree tree, boolean copyAttributes) {

		super(name);

		// get a rooted version of the tree to clone
		FlexibleTree binaryTree = new FlexibleTree(tree, copyAttributes);
		binaryTree.resolveTree();

		// adjust the heights to be compatible with the tip dates and perturb
		// any zero branches.
		MutableTree.Utils.correctHeightsForTips(binaryTree);

		// clone the node structure (this will create the individual parameters)
		Node node = new Node(binaryTree, binaryTree.getRoot());

		internalNodeCount = binaryTree.getInternalNodeCount();
		externalNodeCount = binaryTree.getExternalNodeCount();
		trueInternalNodeCount = internalNodeCount;
		trueExternalNodeCount = externalNodeCount;
		trueNodeCount = internalNodeCount + externalNodeCount;

//		System.err.println(internalNodeCount +"\t"+ externalNodeCount);
		externalNodeCount = initExternalNode;
		while(trueExternalNodeCount > externalNodeCount ){
//			System.err.println(internalNodeCount +"\t"+ externalNodeCount);
			externalNodeCount *= SCALE_TAXA_FACTOR;
		}
		internalNodeCount = externalNodeCount-1;
		nodeCount = internalNodeCount + externalNodeCount;

//		System.err.println(externalNode/Count +"\t"+ internalNodeCount );

		nodes = new Node[nodeCount];
		storedNodes = new Node[nodeCount];

		int i = 0;
		int j = externalNodeCount;

		root = node;

		do {
			node = (Node) Tree.Utils.postorderSuccessor(this, node);
//			System.err.print(node.number +"\t" + node.getHeight());
			if (node.isExternal()) {
				node.number = i;

				nodes[i] = node;
				storedNodes[i] = new Node();
				storedNodes[i].taxon = node.taxon;
				storedNodes[i].number = i;
				i++;
				
			} else {
				node.number = j;

				nodes[j] = node;
				storedNodes[j] = new Node();
				storedNodes[j].number = j;

				j++;
			}
//			System.err.println("\t"+ node.number +"\t" + node.getHeight());
		} while (node != root);
		
//		System.out.println("i"+i +"\tj"+ j);
		for (; i < externalNodeCount; i++) {
			nodes[i] = createEmptyNode(i);
			storedNodes[i] = createEmptyNode(i);
		}
		for (; j < nodeCount; j++) {
			nodes[j] = createEmptyNode(j);
			storedNodes[j] = createEmptyNode(j);
		}
		
		int trueRootIndex = root.getNumber();
		trueRoot = nodes[trueRootIndex];
		
//		System.out.println(externalNodeCount +"\t"+ internalNodeCount );
//		System.out.println(trueExternalNodeCount +"\t"+ trueInternalNodeCount);
//		System.out.println("start: "+j +"\t"+ trueRootIndex );
		int trueRootIndexPlus1 = trueRootIndex+1;
		int internalNodeCountMinus1 = internalNodeCount-1;
		if(trueRootIndexPlus1 != j){
			root = nodes[--j];
			root.addChild(nodes[trueRootIndex]);
//			System.out.println("connect fake root"+root.toString() +"\t"+ root.getChildCount());
			if(j> trueRootIndexPlus1){
				root.addChild(nodes[--j]);

				while(j>trueRootIndexPlus1){
					int otherChildIndex = j-internalNodeCountMinus1;
					nodes[j].addChild(nodes[otherChildIndex]);
					nodes[j].addChild(nodes[--j]);
				}
				int otherChildIndex = j-internalNodeCountMinus1;
//					System.out.println("last j\t"+j +"\t"+otherChildIndex);
				nodes[j].addChild(nodes[otherChildIndex]);
				nodes[j].addChild(nodes[--otherChildIndex]);
				
			}
			else{
				root.addChild(nodes[externalNodeCount-1]);
			}
				
			
		}
		// must be done here to allow programmatic running of BEAST
		setupHeightBounds();
//		System.out.println(root.getNumber() +"\t"+ trueRoot.getNumber());
//		Node newNode = new Node(this, root);
		
	}

	private Node createEmptyNode(int i) {
		Node n = new Node();
		n.number = i;
		n.heightParameter = new Parameter.Default(0.0);
		n.taxon = new Taxon(Integer.toString(i));
		addVariable(n.heightParameter);
		return n;
	}


	/*
	 * New constructor that does NOT alter tree branch length
	 */

	public SSTreeModel(String name, Tree tree, boolean copyAttributes,
			boolean isCorrectHeight) {

		super(name);

		// get a rooted version of the tree to clone
		FlexibleTree binaryTree = new FlexibleTree(tree, copyAttributes);
		binaryTree.resolveTree();

		// adjust the heights to be compatible with the tip dates and perturb
		// any zero branches.
		if (isCorrectHeight) {
			MutableTree.Utils.correctHeightsForTips(binaryTree);
		}

		// clone the node structure (this will create the individual parameters)
		Node node = new Node(binaryTree, binaryTree.getRoot());

		internalNodeCount = binaryTree.getInternalNodeCount();
		externalNodeCount = binaryTree.getExternalNodeCount();

		nodeCount = internalNodeCount + externalNodeCount;

		nodes = new Node[nodeCount];
		storedNodes = new Node[nodeCount];

		int i = 0;
		int j = externalNodeCount;

		root = node;

		do {
			node = (Node) Tree.Utils.postorderSuccessor(this, node);

			if (node.isExternal()) {
				node.number = i;

				nodes[i] = node;
				storedNodes[i] = new Node();
				storedNodes[i].taxon = node.taxon;
				storedNodes[i].number = i;

				i++;
			} else {
				node.number = j;

				nodes[j] = node;
				storedNodes[j] = new Node();
				storedNodes[j].number = j;

				j++;
			}
			
		} while (node != root);

		// must be done here to allow programmatic running of BEAST
		setupHeightBounds();
	}

	boolean heightBoundsSetup = false;

	public void setupHeightBounds() {

		if (heightBoundsSetup) {
			throw new IllegalArgumentException(
					"Node height bounds set up twice");
		}

		for (int i = 0; i < nodeCount; i++) {
			nodes[i].setupHeightBounds();
		}

		heightBoundsSetup = true;
	}

	/**
	 * Push a tree changed event into the event stack.
	 */
	public void pushTreeChangedEvent() {
		pushTreeChangedEvent(new TreeChangedEvent());
	}

	/**
	 * Push a tree changed event into the event stack.
	 */
	public void pushTreeChangedEvent(NodeRef nodeRef) {
		pushTreeChangedEvent(new TreeChangedEvent((Node) nodeRef));
	}

	/**
	 * Push a tree changed event into the event stack.
	 */
	public void pushTreeChangedEvent(Node node, Parameter parameter, int index) {
		pushTreeChangedEvent(new TreeChangedEvent(node, parameter, index));
	}

	/**
	 * Push a tree changed event into the event stack.
	 */
	public void pushTreeChangedEvent(TreeChangedEvent event) {
		if (inEdit) {
			treeChangedEvents.add(event);
		} else {
			listenerHelper.fireModelChanged(this, event);
		}
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// no submodels so nothing to do
	}

	/**
	 * Called when a parameter changes.
	 */
	@Override
	public void handleVariableChangedEvent(Variable variable, int index,
			Parameter.ChangeType type) {
		final Node node = getNodeOfParameter((Parameter) variable);
		if (type == Parameter.ChangeType.ALL_VALUES_CHANGED) {
			// this signals events where values in all dimensions of a parameter
			// is changed.
			pushTreeChangedEvent(new TreeChangedEvent(node,
					(Parameter) variable,
					TreeChangedEvent.CHANGE_IN_ALL_INTERNAL_NODES));
		} else {
			pushTreeChangedEvent(node, (Parameter) variable, index);
		}
	}

	private final List<TreeChangedEvent> treeChangedEvents = new ArrayList<TreeChangedEvent>();

	public boolean hasRates() {
		return hasRates;
	}

	public class TreeChangedEvent {
		static final int CHANGE_IN_ALL_INTERNAL_NODES = -2;

		final Node node;
		final Parameter parameter;
		final int index;

		public TreeChangedEvent() {
			this(null, null, -1);
		}

		public TreeChangedEvent(Node node) {
			this(node, null, -1);
		}

		public TreeChangedEvent(Node node, Parameter parameter, int index) {
			this.node = node;
			this.parameter = parameter;
			this.index = index;
		}

		public int getIndex() {
			return index;
		}

		public Node getNode() {
			return node;
		}

		public Parameter getParameter() {
			return parameter;
		}

		public boolean isTreeChanged() {
			return parameter == null;
		}

		public boolean isNodeChanged() {
			return node != null;
		}

		public boolean isNodeParameterChanged() {
			return parameter != null;
		}

		public boolean isHeightChanged() {
			return parameter == node.heightParameter;
		}

		public boolean isRateChanged() {
			return parameter == node.rateParameter;
		}

		public boolean isTraitChanged(String name) {
			return parameter == node.traitParameters.get(name);
		}

		public boolean areAllInternalHeightsChanged() {
			if (parameter != null) {
				return parameter == node.heightParameter
						&& index == CHANGE_IN_ALL_INTERNAL_NODES;
			}
			return false;
		}

	}

	// *****************************************************************
	// Interface Tree
	// *****************************************************************

	/**
	 * Return the units that this tree is expressed in.
	 */
	@Override
	public Type getUnits() {
		return units;
	}

	/**
	 * Sets the units that this tree is expressed in.
	 */
	@Override
	public void setUnits(Type units) {
		this.units = units;
	}

	/**
	 * @return a count of the number of nodes (internal + external) in this
	 *         tree.
	 */
	@Override
	public int getNodeCount() {
		return nodeCount;
	}

	@Override
	public boolean hasNodeHeights() {
		return true;
	}

	@Override
	public double getNodeHeight(NodeRef node) {
		return ((Node) node).getHeight();
	}

	public final double getNodeHeightUpper(NodeRef node) {
		return ((Node) node).heightParameter.getBounds().getUpperLimit(0);
	}

	public final double getNodeHeightLower(NodeRef node) {
		return ((Node) node).heightParameter.getBounds().getLowerLimit(0);
	}

	/**
	 * @param node
	 * @return the rate parameter associated with this node.
	 */
	@Override
	public double getNodeRate(NodeRef node) {
		if (!hasRates) {
			return 1.0;
		}
		return ((Node) node).getRate();
	}

	@Override
	public Object getNodeAttribute(NodeRef node, String name) {

		if (name.equals("rate")) {
			return getNodeRate(node);
		}

		return null;
	}

	@Override
	public Iterator getNodeAttributeNames(NodeRef node) {
		return new Iterator() {

			int i = 0;
			String[] attributes = { "rate" };

			@Override
			public boolean hasNext() {
				return i < attributes.length;
			}

			@Override
			public Object next() {
				return attributes[i++];
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException(
						"can't remove from this iterator!");
			}
		};
	}

	public boolean hasNodeTraits() {
		return hasTraits;
	}

	public Map<String, Parameter> getTraitMap(NodeRef node) {
		if (!hasTraits)
			throw new IllegalArgumentException(
					"Trait parameters have not been created");
		return ((Node) node).getTraitMap();
	}

	public double getNodeTrait(NodeRef node, String name) {
		if (!hasTraits)
			throw new IllegalArgumentException(
					"Trait parameters have not been created");
		return ((Node) node).getTrait(name);
	}

	public Parameter getNodeTraitParameter(NodeRef node, String name) {
		if (!hasTraits)
			throw new IllegalArgumentException(
					"Trait parameters have not been created");
		return ((Node) node).getTraitParameter(name);
	}

	public double[] getMultivariateNodeTrait(NodeRef node, String name) {
		if (!hasTraits)
			throw new IllegalArgumentException(
					"Trait parameters have not been created");
		return ((Node) node).getMultivariateTrait(name);
	}

	public final void swapAllTraits(NodeRef node1, NodeRef node2) {
		if (!hasTraits)
			throw new IllegalArgumentException(
					"Trait parameters have not been created");
		swapAllTraits((Node) node1, (Node) node2);
	}

	@Override
	public Taxon getNodeTaxon(NodeRef node) {
		return ((Node) node).taxon;
	}

	public void setNodeTaxon(NodeRef node, Taxon taxon) {
		((Node) node).taxon = taxon;
	}

	@Override
	public boolean isExternal(NodeRef node) {
		return ((Node) node).isExternal();
	}

	@Override
	public boolean isRoot(NodeRef node) {
		return (node == root);
	}

	@Override
	public int getChildCount(NodeRef node) {
		return ((Node) node).getChildCount();
	}

	@Override
	public NodeRef getChild(NodeRef node, int i) {
		return ((Node) node).getChild(i);
	}

	@Override
	public NodeRef getParent(NodeRef node) {
		return ((Node) node).parent;
	}

	@Override
	public boolean hasBranchLengths() {
		return true;
	}

	@Override
	public double getBranchLength(NodeRef node) {
		NodeRef parent = getParent(node);
		if (parent == null) {
			return 0.0;
		}

		return getNodeHeight(parent) - getNodeHeight(node);
	}

	@Override
	public NodeRef getExternalNode(int i) {
		return nodes[i];
	}

	@Override
	public NodeRef getInternalNode(int i) {
		return nodes[i + externalNodeCount];
	}

	@Override
	public NodeRef getNode(int i) {
		return nodes[i];
	}

	public NodeRef[] getNodes() {
		return nodes;
	}

	/**
	 * Returns the number of external nodes.
	 */
	@Override
	public int getExternalNodeCount() {
		return externalNodeCount;
	}

	/**
	 * Returns the ith internal node.
	 */
	@Override
	public int getInternalNodeCount() {
		return internalNodeCount;
	}

	/**
	 * Returns the root node of this tree.
	 */
	@Override
	public NodeRef getRoot() {
		return root;
	}

	// *****************************************************************
	// Interface MutableTree
	// *****************************************************************

	/**
	 * Set a new node as root node.
	 */
	@Override
	public final void setRoot(NodeRef newRoot) {

		if (!inEdit)
			throw new RuntimeException(
					"Must be in edit transaction to call this method!");

		root = (Node) newRoot;

		// We shouldn't need this because the addChild will already have fired
		// appropriate events.
		pushTreeChangedEvent(root);
	}

	@Override
	public void addChild(NodeRef p, NodeRef c) {

		if (!inEdit)
			throw new RuntimeException(
					"Must be in edit transaction to call this method!");

		Node parent = (Node) p;
		Node child = (Node) c;
		if (parent.hasChild(child))
			throw new IllegalArgumentException("Child already exists in parent");

		parent.addChild(child);
		pushTreeChangedEvent(parent);
	}

	@Override
	public void removeChild(NodeRef p, NodeRef c) {

		if (!inEdit)
			throw new RuntimeException(
					"Must be in edit transaction to call this method!");

		Node parent = (Node) p;
		Node child = (Node) c;

		parent.removeChild(child);
	}

	@Override
	public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {
		throw new RuntimeException("Unimplemented");
	}

	private Node oldRoot;

	@Override
	public boolean beginTreeEdit() {
		if (inEdit)
			throw new RuntimeException("Alreading in edit transaction mode!");

		oldRoot = root;

		inEdit = true;

		return false;
	}

	@Override
	public void endTreeEdit() {
		if (!inEdit)
			throw new RuntimeException("Not in edit transaction mode!");

		inEdit = false;

		if (root != oldRoot) {
			swapParameterObjects(oldRoot, root);
		}

		if (TEST_NODE_BOUNDS) {
			try {
				checkTreeIsValid();
			} catch (InvalidTreeException ite) {
				throw new RuntimeException(ite.getMessage());
			}
		}

		for (TreeChangedEvent treeChangedEvent : treeChangedEvents) {
			listenerHelper.fireModelChanged(this, treeChangedEvent);
		}
		treeChangedEvents.clear();
	}

	public void checkTreeIsValid() throws MutableTree.InvalidTreeException {
		for (Node node : nodes) {
			if (!node.heightParameter.isWithinBounds()) {
				throw new InvalidTreeException("height parameter out of bounds");
			}
		}
	}

	@Override
	public void setNodeHeight(NodeRef n, double height) {
		((Node) n).setHeight(height);
	}

	@Override
	public void setNodeRate(NodeRef n, double rate) {
		if (!hasRates)
			throw new IllegalArgumentException(
					"Rate parameters have not been created");
		((Node) n).setRate(rate);

	}

	public void setNodeTrait(NodeRef n, String name, double value) {
		if (!hasTraits)
			throw new IllegalArgumentException(
					"Trait parameters have not been created");
		((Node) n).setTrait(name, value);
	}

	public void setMultivariateTrait(NodeRef n, String name, double[] value) {
		if (!hasTraits)
			throw new IllegalArgumentException(
					"Trait parameters have not been created");
		((Node) n).setMultivariateTrait(name, value);
	}

	@Override
	public void setBranchLength(NodeRef node, double length) {
		throw new UnsupportedOperationException(
				"TreeModel cannot have branch lengths set");
	}

	@Override
	public void setNodeAttribute(NodeRef node, String name, Object value) {
		throw new UnsupportedOperationException(
				"TreeModel does not use NodeAttributes");
	}

	// *****************************************************************
	// Interface ModelComponent
	// *****************************************************************

	/**
	 * Store current state
	 */
	@Override
	protected void storeState() {

		copyNodeStructure(storedNodes);
		storedRootNumber = root.getNumber();

	}

	/**
	 * Restore the stored state
	 */
	@Override
	protected void restoreState() {

		Node[] tmp = storedNodes;
		storedNodes = nodes;
		nodes = tmp;

		root = nodes[storedRootNumber];
	}

	/**
	 * accept the stored state
	 */
	@Override
	protected void acceptState() {
	} // nothing to do

	/**
	 * Copies the node connections from this TreeModel's nodes array to the
	 * destination array. Basically it connects up the nodes in destination in
	 * the same way as this TreeModel is set up. This method is package private.
	 */
	void copyNodeStructure(Node[] destination) {

		if (nodes.length != destination.length) {
			throw new IllegalArgumentException(
					"Node arrays are of different lengths");
		}

		for (int i = 0, n = nodes.length; i < n; i++) {
			Node node0 = nodes[i];
			Node node1 = destination[i];

			// the parameter values are automatically stored and restored
			// just need to keep the links
			node1.heightParameter = node0.heightParameter;
			node1.rateParameter = node0.rateParameter;
			node1.traitParameters = node0.traitParameters;

			if (node0.parent != null) {
				node1.parent = storedNodes[node0.parent.getNumber()];
			} else {
				node1.parent = null;
			}

			if (node0.leftChild != null) {
				node1.leftChild = storedNodes[node0.leftChild.getNumber()];
			} else {
				node1.leftChild = null;
			}

			if (node0.rightChild != null) {
				node1.rightChild = storedNodes[node0.rightChild.getNumber()];
			} else {
				node1.rightChild = null;
			}
		}
	}

	/**
	 * @return the number of statistics of this component.
	 */
	@Override
	public int getStatisticCount() {
		return super.getStatisticCount() + 1;
	}

	/**
	 * @return the ith statistic of the component
	 */
	@Override
	public Statistic getStatistic(int i) {
		if (i == super.getStatisticCount())
			return root.heightParameter;
		return super.getStatistic(i);
	}

	// public String getModelComponentName() {
	// return TREE_MODEL;
	// }

	// **************************************************************
	// TaxonList IMPLEMENTATION
	// **************************************************************

	/**
	 * @return a count of the number of taxa in the list.
	 */
	@Override
	public int getTaxonCount() {
		return getExternalNodeCount();
	}

	/**
	 * @return the ith taxon in the list.
	 */
	@Override
	public Taxon getTaxon(int taxonIndex) {
		return ((Node) getExternalNode(taxonIndex)).taxon;
	}

	/**
	 * @return the ID of the taxon of the ith external node. If it doesn't have
	 *         a taxon, returns the ID of the node itself.
	 */
	@Override
	public String getTaxonId(int taxonIndex) {
		Taxon taxon = getTaxon(taxonIndex);
		if (taxon != null) {
			return taxon.getId();
		} else {
			return null;
		}
	}

	/**
	 * returns the index of the taxon with the given id.
	 */
	@Override
	public int getTaxonIndex(String id) {
		for (int i = 0, n = getTaxonCount(); i < n; i++) {
			if (getTaxonId(i).equals(id))
				return i;
		}
		return -1;
	}

	/**
	 * returns the index of the given taxon.
	 */
	@Override
	public int getTaxonIndex(Taxon taxon) {
		for (int i = 0, n = getTaxonCount(); i < n; i++) {
			if (getTaxon(i) == taxon)
				return i;
		}
		return -1;
	}

	@Override
	public List<Taxon> asList() {
		List<Taxon> taxa = new ArrayList<Taxon>();
		for (int i = 0, n = getTaxonCount(); i < n; i++) {
			taxa.add(getTaxon(i));
		}
		return taxa;
	}

	@Override
	public Iterator<Taxon> iterator() {
		return new Iterator<Taxon>() {
			private int index = -1;

			@Override
			public boolean hasNext() {
				return index < getTaxonCount() - 1;
			}

			@Override
			public Taxon next() {
				index++;
				return getTaxon(index);
			}

			@Override
			public void remove() { /* do nothing */
			}
		};
	}

	/**
	 * @param taxonIndex
	 *            the index of the taxon whose attribute is being fetched.
	 * @param name
	 *            the name of the attribute of interest.
	 * @return an object representing the named attributed for the taxon of the
	 *         given external node. If the node doesn't have a taxon then the
	 *         nodes own attribute is returned.
	 */
	@Override
	public Object getTaxonAttribute(int taxonIndex, String name) {
		Taxon taxon = getTaxon(taxonIndex);
		if (taxon != null) {
			return taxon.getAttribute(name);
		}
		return null;
	}

	// **************************************************************
	// MutableTaxonList IMPLEMENTATION
	// **************************************************************

	@Override
	public int addTaxon(Taxon taxon) {
		throw new IllegalArgumentException("Cannot add taxon to a TreeModel");
	}

	@Override
	public boolean removeTaxon(Taxon taxon) {
		throw new IllegalArgumentException("Cannot add taxon to a TreeModel");
	}

	@Override
	public void setTaxonId(int taxonIndex, String id) {
		throw new IllegalArgumentException("Cannot set taxon id in a TreeModel");
	}

	@Override
	public void setTaxonAttribute(int taxonIndex, String name, Object value) {
		throw new IllegalArgumentException(
				"Cannot set taxon attribute in a TreeModel");
	}

	@Override
	public void addMutableTreeListener(MutableTreeListener listener) {
	} // Do nothing at the moment

	@Override
	public void addMutableTaxonListListener(MutableTaxonListListener listener) {
	} // Do nothing at the moment

	// **************************************************************
	// Identifiable IMPLEMENTATION
	// **************************************************************

	private String id = null;

	/**
	 * @return the id.
	 */
	@Override
	public String getId() {
		return id;
	}

	/**
	 * Sets the id.
	 */
	@Override
	public void setId(String id) {
		this.id = id;
	}

	// **************************************************************
	// Attributable IMPLEMENTATION
	// **************************************************************

	private Attributable.AttributeHelper treeAttributes = null;

	/**
	 * Sets an named attribute for this object.
	 * 
	 * @param name
	 *            the name of the attribute.
	 * @param value
	 *            the new value of the attribute.
	 */
	@Override
	public void setAttribute(String name, Object value) {
		if (treeAttributes == null)
			treeAttributes = new Attributable.AttributeHelper();
		treeAttributes.setAttribute(name, value);
	}

	/**
	 * @param name
	 *            the name of the attribute of interest.
	 * @return an object representing the named attributed for this object.
	 */
	@Override
	public Object getAttribute(String name) {
		if (treeAttributes == null)
			return null;
		else
			return treeAttributes.getAttribute(name);
	}

	/**
	 * @return an iterator of the attributes that this object has.
	 */
	@Override
	public Iterator<String> getAttributeNames() {
		if (treeAttributes == null)
			return null;
		else
			return treeAttributes.getAttributeNames();
	}

	/**
	 * @return a string containing a newick representation of the tree
	 */
	public final String getNewick() {
		return Tree.Utils.newick(this);
	}

	/**
	 * @return a string containing a newick representation of the tree
	 */
	@Override
	public String toString() {
		return getNewick();
	}

	@Override
	public Tree getCopy() {
		throw new UnsupportedOperationException(
				"please don't call this function");
	}

	// **************************************************************
	// XMLElement IMPLEMENTATION
	// **************************************************************

	@Override
	public Element createElement(Document document) {
		throw new RuntimeException("Not implemented yet");
	}

	// ***********************************************************************
	// Private methods
	// ***********************************************************************

	/**
	 * @return the node that this parameter is a member of
	 */
	public Node getNodeOfParameter(Parameter parameter) {

		if (parameter == null)
			throw new IllegalArgumentException("Parameter is null!");

		for (Node node : nodes) {
			if (node.heightParameter == parameter) {
				return node;
			}
		}

		if (hasRates) {
			for (Node node : nodes) {
				if (node.rateParameter == parameter) {
					return node;
				}
			}
		}
		if (hasTraits) {
			for (Node node : nodes) {
				if (node.traitParameters.containsValue(parameter)) {
					return node;
				}
			}
		}
		throw new RuntimeException("Parameter not found in any nodes:"
				+ parameter.getId() + " " + parameter.hashCode());
		// assume it is a trait parameter and return null
		// return null;
	}

	/**
	 * Get the root height parameter. Is private because it can only be called
	 * by the XMLParser
	 */
	public Parameter getRootHeightParameter() {

		return root.heightParameter;
	}

	/**
	 * @return the relevant node height parameter. Is private because it can
	 *         only be called by the XMLParser
	 */
	public Parameter createNodeHeightsParameter(boolean rootNode,
			boolean internalNodes, boolean leafNodes) {

		if (!rootNode && !internalNodes && !leafNodes) {
			throw new IllegalArgumentException(
					"At least one of rootNode, internalNodes or leafNodes must be true");
		}

		CompoundParameter parameter = new CompoundParameter("nodeHeights("
				+ getId() + ")");

		for (int i = externalNodeCount; i < nodeCount; i++) {
			if ((rootNode && nodes[i] == root)
					|| (internalNodes && nodes[i] != root)) {
				parameter.addParameter(nodes[i].heightParameter);
			}
		}

		if (leafNodes) {
			for (int i = 0; i < externalNodeCount; i++) {
				parameter.addParameter(nodes[i].heightParameter);
			}
		}

		return parameter;
	}

	public Parameter getLeafHeightParameter(NodeRef node) {

		if (!isExternal(node)) {
			throw new RuntimeException(
					"only root and leaves can be used with setNodeHeightParameter");
		}

		return nodes[node.getNumber()].heightParameter;
	}

	/**
	 * @return the relevant node rate parameter. Is private because it can only
	 *         be called by the XMLParser
	 */
	public Parameter createNodeRatesParameter(double[] initialValues,
			boolean rootNode, boolean internalNodes, boolean leafNodes) {

		if (!rootNode && !internalNodes && !leafNodes) {
			throw new IllegalArgumentException(
					"At least one of rootNode, internalNodes or leafNodes must be true");
		}

		CompoundParameter parameter = new CompoundParameter("nodeRates("
				+ getId() + ")");

		hasRates = true;

		for (int i = externalNodeCount; i < nodeCount; i++) {
			nodes[i].createRateParameter(initialValues);
			if ((rootNode && nodes[i] == root)
					|| (internalNodes && nodes[i] != root)) {
				parameter.addParameter(nodes[i].rateParameter);
			}
		}

		for (int i = 0; i < externalNodeCount; i++) {
			nodes[i].createRateParameter(initialValues);
			if (leafNodes) {
				parameter.addParameter(nodes[i].rateParameter);
			}
		}

		return parameter;
	}

	public Parameter createNodeTraitsParameter(String name,
			double[] initialValues) {
		return createNodeTraitsParameter(name, initialValues.length,
				initialValues, true, true, true, true);
	}

	/**
	 * Create a node traits parameter. Is private because it can only be called
	 * by the XMLParser
	 */
	public Parameter createNodeTraitsParameter(String name, int dim,
			double[] initialValues, boolean rootNode, boolean internalNodes,
			boolean leafNodes, boolean firesTreeEvents) {

		if (!rootNode && !internalNodes && !leafNodes) {
			throw new IllegalArgumentException(
					"At least one of rootNode, internalNodes or leafNodes must be true");
		}

		CompoundParameter parameter = new CompoundParameter(name);

		hasTraits = true;

		for (int i = externalNodeCount; i < nodeCount; i++) {
			nodes[i].createTraitParameter(name, dim, initialValues,
					firesTreeEvents);
			if ((rootNode && nodes[i] == root)
					|| (internalNodes && nodes[i] != root)) {
				parameter.addParameter(nodes[i].getTraitParameter(name));
			}
		}

		for (int i = 0; i < externalNodeCount; i++) {
			nodes[i].createTraitParameter(name, dim, initialValues,
					firesTreeEvents);
			if (leafNodes) {
				parameter.addParameter(nodes[i].getTraitParameter(name));
			}
		}

		return parameter;
	}

	private void swapAllTraits(Node n1, Node n2) {

		for (Map.Entry<String, Parameter> entry : n1.traitParameters.entrySet()) {
			Parameter p1 = n1.traitParameters.get(entry.getKey());
			Parameter p2 = n2.traitParameters.get(entry.getKey());
			final int dim = p1.getDimension();
			for (int i = 0; i < dim; i++) {
				double transfer = p1.getParameterValue(i);
				p1.setParameterValue(i, p2.getParameterValue(i));
				p2.setParameterValue(i, transfer);
			}

		}

	}

	/**
	 * This method swaps the parameter objects of the two nodes but maintains
	 * the values in each node. This method is used to ensure that root node of
	 * the tree always has the same parameter object.
	 */
	private void swapParameterObjects(Node n1, Node n2) {

		double height1 = n1.getHeight();
		double height2 = n2.getHeight();

		double rate1 = 1.0, rate2 = 1.0;

		if (hasRates) {
			rate1 = n1.getRate();
			rate2 = n2.getRate();
		}

		// swap all trait parameters

		if (hasTraits) {
			Map<String, Parameter> traits1 = new HashMap<String, Parameter>();
			Map<String, Parameter> traits2 = new HashMap<String, Parameter>();

			traits1.putAll(n1.traitParameters);
			traits2.putAll(n2.traitParameters);

			Map<String, Parameter> temp = n1.traitParameters;
			n1.traitParameters = n2.traitParameters;
			n2.traitParameters = temp;

			for (Map.Entry<String, Parameter> entry : traits1.entrySet()) {
				n1.traitParameters.get(entry.getKey())
						.setParameterValueQuietly(0,
								entry.getValue().getParameterValue(0));
			}
			for (Map.Entry<String, Parameter> entry : traits2.entrySet()) {
				n2.traitParameters.get(entry.getKey())
						.setParameterValueQuietly(0,
								entry.getValue().getParameterValue(0));
			}
		}

		Parameter temp = n1.heightParameter;
		n1.heightParameter = n2.heightParameter;
		n2.heightParameter = temp;

		if (hasRates) {
			temp = n1.rateParameter;
			n1.rateParameter = n2.rateParameter;
			n2.rateParameter = temp;
		}

		n1.heightParameter.setParameterValueQuietly(0, height1);
		n2.heightParameter.setParameterValueQuietly(0, height2);

		if (hasRates) {
			n1.rateParameter.setParameterValueQuietly(0, rate1);
			n2.rateParameter.setParameterValueQuietly(0, rate2);
		}
	}

	// **************************************************************
	// Private inner classes
	// **************************************************************

	public class Node implements NodeRef {

		public Node parent;
		public Node leftChild, rightChild;
		private int number;
		public Parameter heightParameter;
		public Parameter rateParameter = null;
		// public Parameter traitParameter = null;
		public Taxon taxon = null;

		Map<String, Parameter> traitParameters = new HashMap<String, Parameter>();

		public Node() {
			parent = null;
			leftChild = rightChild = null;
			heightParameter = null;
			number = 0;
			taxon = null;
		}

		/**
		 * constructor used to clone a node and all children
		 */
		public Node(Tree tree, NodeRef node) {
			parent = null;
			leftChild = rightChild = null;

			heightParameter = new Parameter.Default(tree.getNodeHeight(node));
			addVariable(heightParameter);

			number = node.getNumber();
			taxon = tree.getNodeTaxon(node);
			heightParameter.setId("" + number);
			for (int i = 0; i < tree.getChildCount(node); i++) {
				addChild(new Node(tree, tree.getChild(node, i)));
			}
		}

		public final void setupHeightBounds() {
			heightParameter.addBounds(new NodeHeightBounds(heightParameter));
		}

		public final void createRateParameter(double[] initialValues) {
			if (rateParameter == null) {
				if (initialValues != null) {
					rateParameter = new Parameter.Default(initialValues[0]);
				} else {
					rateParameter = new Parameter.Default(1.0);
				}
				if (isRoot()) {
					rateParameter.setId("root.rate");
				} else if (isExternal()) {
					rateParameter.setId(getTaxonId(getNumber()) + ".rate");
				} else {
					rateParameter.setId("node" + getNumber() + ".rate");
				}
				rateParameter.addBounds(new Parameter.DefaultBounds(
						Double.POSITIVE_INFINITY, 0.0, 1));
				addVariable(rateParameter);
			}
		}

		public final void createTraitParameter(String name,
				double[] initialValues, boolean firesTreeEvents) {
			createTraitParameter(name, initialValues.length, initialValues,
					firesTreeEvents);
		}

		public final void createTraitParameter(String name, int dim,
				double[] initialValues, boolean firesTreeEvents) {

			if (!traitParameters.containsKey(name)) {

				Parameter trait = new Parameter.Default(dim);
				if (isRoot()) {
					trait.setId("root." + name);
				} else if (isExternal()) {
					trait.setId(getTaxonId(getNumber()) + "." + name);
				} else {
					trait.setId("node" + getNumber() + "." + name);
				}
				trait.addBounds(new Parameter.DefaultBounds(
						Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, dim));

				if (initialValues != null && initialValues.length > 0) {
					for (int i = 0; i < dim; i++) {
						if (initialValues.length == dim) {
							trait.setParameterValue(i, initialValues[i]);
						} else {
							trait.setParameterValue(i, initialValues[0]);
						}
					}
				}

				traitParameters.put(name, trait);

				if (firesTreeEvents) {
					addVariable(trait);
				}
			}
		}

		public final double getHeight() {
			return heightParameter.getParameterValue(0);
		}

		public final double getRate() {
			return rateParameter.getParameterValue(0);
		}

		public final double getTrait(String name) {
			return traitParameters.get(name).getParameterValue(0);
		}

		public final double[] getMultivariateTrait(String name) {
			return traitParameters.get(name).getParameterValues();
		}

		public final Map<String, Parameter> getTraitMap() {
			return traitParameters;
		}

		public final void setHeight(double height) {
			heightParameter.setParameterValue(0, height);
		}

		public final void setRate(double rate) {
			// System.out.println("Rate set for parameter " +
			// rateParameter.getParameterName());
			rateParameter.setParameterValue(0, rate);
		}

		public final void setTrait(String name, double trait) {
			// System.out.println("Trait set for parameter " +
			// traitParameter.getParameterName());
			traitParameters.get(name).setParameterValue(0, trait);
		}

		public final void setMultivariateTrait(String name, double[] trait) {
			int dim = trait.length;
			for (int i = 0; i < dim; i++)
				traitParameters.get(name).setParameterValue(i, trait[i]);
		}

		@Override
		public int getNumber() {
			return number;
		}

		@Override
		public void setNumber(int n) {
			number = n;
		}

		/**
		 * Returns the number of children this node has.
		 */
		public final int getChildCount() {
			int n = 0;
			if (leftChild != null)
				n++;
			if (rightChild != null)
				n++;
			return n;
		}

		public Node getChild(int n) {
			if (n == 0)
				return leftChild;
			if (n == 1)
				return rightChild;
			throw new IllegalArgumentException(
					"TreeModel.Nodes can only have 2 children");
		}

		public boolean hasChild(Node node) {
			return (leftChild == node || rightChild == node);
		}

		/**
		 * add new child node
		 * 
		 * @param node
		 *            new child node
		 */
		public void addChild(Node node) {
			if (leftChild == null) {
				leftChild = node;
			} else if (rightChild == null) {
				rightChild = node;
			} else {
				throw new IllegalArgumentException(
						"TreeModel.Nodes can only have 2 children");
			}
			node.parent = this;
		}

		/**
		 * remove child
		 * 
		 * @param node
		 *            child to be removed
		 */
		public Node removeChild(Node node) {
			if (leftChild == node) {
				leftChild = null;
			} else if (rightChild == node) {
				rightChild = null;
			} else {
				throw new IllegalArgumentException("Unknown child node");
			}
			node.parent = null;
			return node;
		}

		/**
		 * remove child
		 * 
		 * @param n
		 *            number of child to be removed
		 */
		public Node removeChild(int n) {
			Node node;
			if (n == 0) {
				node = leftChild;
				leftChild = null;
			} else if (n == 1) {
				node = rightChild;
				rightChild = null;
			} else {
				throw new IllegalArgumentException(
						"TreeModel.Nodes can only have 2 children");
			}
			node.parent = null;
			return node;
		}

		public boolean hasNoChildren() {
			return (leftChild == null && rightChild == null);
		}

		public boolean isExternal() {
			return hasNoChildren();
		}

		public boolean isRoot() {
			return (parent == null);
		}

		@Override
		public String toString() {
			return "node " + number + ", height=" + getHeight()
					+ (taxon != null ? ": " + taxon.getId() : "");
		}

		public Parameter getTraitParameter(String name) {
			return traitParameters.get(name);
		}
	}

	/**
	 * This class provides bounds for parameters that represent a node height in
	 * this tree model.
	 */
	private class NodeHeightBounds implements Bounds<Double> {

		public NodeHeightBounds(Parameter parameter) {
			nodeHeightParameter = parameter;
		}

		@Override
		public Double getUpperLimit(int i) {

			Node node = getNodeOfParameter(nodeHeightParameter);
			if (node.isRoot()) {
				return Double.POSITIVE_INFINITY;
			} else {
				return node.parent.getHeight();
			}
		}

		@Override
		public Double getLowerLimit(int i) {

			Node node = getNodeOfParameter(nodeHeightParameter);
			if (node.isExternal()) {
				return 0.0;
			} else {
				return Math.max(node.leftChild.getHeight(),
						node.rightChild.getHeight());
			}
		}

		@Override
		public int getBoundsDimension() {
			return 1;
		}

		private Parameter nodeHeightParameter = null;
	}

	// ***********************************************************************
	// Private members
	// ***********************************************************************

	/**
	 * root node
	 */
	private Node root = null;
	private int storedRootNumber;

	/**
	 * list of internal nodes (including root)
	 */
	private Node[] nodes = null;
	private Node[] storedNodes = null;

	/**
	 * number of nodes (including root and tips)
	 */
	private final int nodeCount;

	/**
	 * number of external nodes
	 */
	private int externalNodeCount;

	/**
	 * number of internal nodes (including root)
	 */
	private int internalNodeCount;

	/**
	 * holds the units of the trees branches.
	 */
	private Type units = Type.SUBSTITUTIONS;

	private boolean inEdit = false;

	private boolean hasRates = false;
	private boolean hasTraits = false;

	public SSTreeModel insertNode() {

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
		newExternalNode.number = trueExternalNodeCount;
		newExternalNode.heightParameter = new Parameter.Default(newExternalHeight);
		newExternalNode.taxon = new Taxon("XXXhap_"+trueExternalNodeCount);
		addVariable(newExternalNode.heightParameter);
//		newExternalNode.setHeight(newExternalHeight);
		nodes[trueExternalNodeCount] = newExternalNode;
		
//		Node newInternalNode = createEmptyNode(trueNodeCount+1);
		Node newInternalNode = new Node();
		newInternalNode.number = trueNodeCount+1;
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
		SSTreeModel a = new SSTreeModel(st);
		
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
