package srp.dbgraph;

import java.util.ArrayList;
import java.util.Iterator;

import dr.evolution.util.Taxon;
import dr.util.Identifiable;

public class CompatibleSets implements Iterable<CompatibleNode>  {

	private ArrayList<CompatibleNode> allSets;
	private int maxNodeIndex;
	
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

	protected int getNodeCount() {
		return allSets.size();
	}

	protected CompatibleNode getCompatibleNode(int index) {
		return allSets.get(index);
	}
}