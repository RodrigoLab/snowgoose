package srp.dbgraph;

import java.util.Hashtable;

public class DeBruijnGraph {
	
	private Hashtable<Integer, String> allNodes;
	private Hashtable<Integer, Integer> allEdges;
	private int maxNode;
	public DeBruijnGraph() {
		allNodes = new Hashtable<Integer, String>();
		allEdges = new Hashtable<Integer, Integer>();
		maxNode = 0;
	}
	
	public void addNode(int node, String seqs) {
		allNodes.put(node, seqs);
		if(node> maxNode){
			maxNode = node;
		}
	}

	public void addEdge(int node, int node2) {
		allEdges.put(node, node2);
		
	}
	
	public void test(){
		for (int i = 0; i < maxNode; i++) {
			System.out.println(i +"\t"+ allEdges.get(i) +"\t"+ allNodes.get(i));
			
		}
	}

}
