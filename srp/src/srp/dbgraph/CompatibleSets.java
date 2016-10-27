package srp.dbgraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

public class CompatibleSets implements Iterable<CompatibleNode>  {

	private ArrayList<CompatibleNode> allSets;
	private CompatibleNode[] allCNodesArray;
	private int maxNodeIndex;
	private int nodeCount;
	
	public CompatibleSets() {
		
		allSets = new ArrayList<CompatibleNode>();
	}

	public void addCompatibleNode(CompatibleNode cNode){
		allSets.add(cNode);
		int nodeIndex = cNode.getNodeIndex();
		if(nodeIndex> maxNodeIndex){
			maxNodeIndex = nodeIndex;
		}
	}
	
	public int getMaxNodeIndex(){
		return maxNodeIndex;
	}
	
	@Override
	@Deprecated
	public Iterator<CompatibleNode> iterator() {
	    return new Iterator<CompatibleNode>() {
	        private int index = -1;
	
	        @Override
			public boolean hasNext() {
	            return index < getNodeCount() - 1;
	        }
	

			@Override
			public CompatibleNode next() {
	            index++;
	            return getCompatibleNode(index);
	        }
	
	        @Override
			public void remove() { /* do nothing */ }
	    };
	}

	public int getNodeCount() {
		return nodeCount;
	}

//	protected int getNodeCount() {
//		return allSets.size();
//	}

	public CompatibleNode getCompatibleNode(int index) {
		return allCNodesArray[index];
	}

	public void preprocess() {
		for (CompatibleNode cNode : allSets) {
			cNode.proprocess();
		}
		nodeCount = allSets.size();
		allCNodesArray = new CompatibleNode[nodeCount];
		for (int i = 0; i < allCNodesArray.length; i++) {
			CompatibleNode compatibleNode = allSets.get(i);
			int arrayIndex = compatibleNode.getNodeIndex();
			allCNodesArray[arrayIndex] =compatibleNode;  
		}
		
	}
}