package srp.haplotypes;

import java.util.ArrayList;
import java.util.List;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.tree.MutableTree.InvalidTreeException;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeModel.Node;
import dr.evomodel.tree.TreeModel.TreeChangedEvent;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.inference.model.VariableListener;
import dr.inference.model.Model.ListenerHelper;
@Deprecated
public class AlignmentModel extends AbstractModel{
	
	public static final String ALIGNMENT_MODEL = "alignmentModel";
	
	private List<String> alignment = new ArrayList<String>();
	private List<String> storedAlignment = new ArrayList<String>();
	
	private boolean inEdit;
	
	public AlignmentModel(String name) {
		super(name);
	}
	
	public AlignmentModel(ArrayList<String> alignment) {
		this(ALIGNMENT_MODEL);
		this.alignment = alignment;
		
	}
	public AlignmentModel(Alignment alignment) {
		this(ALIGNMENT_MODEL);
		
		ArrayList<String> haplotypes = new ArrayList<String>();
		for (int i = 0; i < alignment.getSequenceCount(); i++) {
			haplotypes.add(alignment.getAlignedSequenceString(i));
		}
		
		this.alignment = haplotypes;
		
	}
	

	

	
//	public void endTreeEdit() { //TreeModel
//        if (!inEdit) throw new RuntimeException("Not in edit transaction mode!");
//
//        inEdit = false;
//
//        if (root != oldRoot) {
//            swapParameterObjects(oldRoot, root);
//        }
//
//        if (TEST_NODE_BOUNDS) {
//            try {
//                checkTreeIsValid();
//            } catch (InvalidTreeException ite) {
//                throw new RuntimeException(ite.getMessage());
//            }
//        }
//
//        for (TreeChangedEvent treeChangedEvent : treeChangedEvents) {
//            ListenerHelper.fireModelChanged(this, treeChangedEvent);
//        }
//        treeChangedEvents.clear();
//    }
    
 
	
//    public void setParameterValue(int i, double val) {//Parameter.Default
//        values[i] = val;
//        fireParameterChangedEvent(i, Parameter.ChangeType.VALUE_CHANGED);
//    }
//    public void fireParameterChangedEvent(int index, Parameter.ChangeType type) {
//        if (listeners != null) {
//            for (VariableListener listener : listeners) {
//                listener.variableChangedEvent(this, index, type);
//            }
//        }
//    }

/*

	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// no submodels so nothing to do
	}


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
*/
	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// no submodels so nothing to do
		// copied from TreeModel!!
	}



	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		// TODO Auto-generated method stub
		
	}



	@Override
	protected void storeState() {
		// TODO Auto-generated method stub

//		copyNodeStructure(storedNodes);
//		storedRootNumber = root.getNumber();
		storedAlignment.clear();
		for (String s : alignment) {
			storedAlignment.add(s);
		}
	}



	@Override
	protected void restoreState() {
		// TODO Auto-generated method stub

		List<String> tmp = storedAlignment;
		storedAlignment = alignment;
		alignment = tmp;
//		Node[] tmp = storedNodes;
//		storedNodes = nodes;
//		nodes = tmp;
//
//		root = nodes[storedRootNumber];
		
		
	}



	@Override
	protected void acceptState() {
		// TODO Auto-generated method stub
		 // nothing to do
	}
	
    public boolean beginAlignmentEdit() {
        if (inEdit) throw new RuntimeException("Alreading in edit transaction mode!");

//        oldRoot = root;

        inEdit = true;

        return false;
    }

    public void endAlignmentEdit() {
        if (!inEdit) throw new RuntimeException("Not in edit transaction mode!");

        inEdit = false;

//        if (root != oldRoot) {
//            swapParameterObjects(oldRoot, root);
//        }

//        if (TEST_NODE_BOUNDS) {
//            try {
//                checkTreeIsValid();
//            } catch (InvalidTreeException ite) {
//                throw new RuntimeException(ite.getMessage());
//            }
//        }

        for (AlignmentChangedEvent alignmentChangedEvent : alignmentChangedEvents) {
            listenerHelper.fireModelChanged(this, alignmentChangedEvent);
        }
        alignmentChangedEvents.clear();
    }
    private final List<AlignmentChangedEvent> alignmentChangedEvents = new ArrayList<AlignmentChangedEvent>();


    public class AlignmentChangedEvent {
        static final int CHANGE_IN_ALL_INTERNAL_NODES = -2;

        final Node node;
        final Parameter parameter;
        final int index;

        public AlignmentChangedEvent() {
            this(null, null, -1);
        }

//        public AlignmentChangedEvent(Node node) {
//            this(node, null, -1);
//        }

        public AlignmentChangedEvent(Node node, Parameter parameter, int index) {
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


        public boolean areAllInternalHeightsChanged() {
            if (parameter != null) {
                return parameter == node.heightParameter && index == CHANGE_IN_ALL_INTERNAL_NODES;
            }
            return false;
        }

    }
}
