package srp.dbgraph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;

import srp.haplotypes.Haplotype;
import srp.haplotypes.HaplotypeModel;

public class DeBruijnGraph {
//	private OpenIntIntHashMap
	private HashMap<Integer, String> allNodes;
	private HashMap<Integer, Integer> allLength;
	private HashMap<Integer, Integer> allEdges;
//	private int maxNode;
//	private CompatibleSets allSets;
	
	public DeBruijnGraph() {
		allNodes = new HashMap<Integer, String>();
		allEdges = new HashMap<Integer, Integer>();
//		maxNode = 0;
	}
	
	public HashMap<Integer, String> getAllNodes() {
		return allNodes;
	}

	public HashMap<Integer, Integer> getAllLength() {
		return allLength;
	}

	public HashMap<Integer, Integer> getAllEdges() {
		return allEdges;
	}

	public void addNode(int node, String seqs) {
		allNodes.put(node, seqs);
		
		allLength.put(node, seqs.length());
	}

	public void addEdge(int node, int node2) {
		allEdges.put(node, node2);
		
	}
	
	public void postPrecossDeBruijnGraph(){
		
		
	}
	
	
}
