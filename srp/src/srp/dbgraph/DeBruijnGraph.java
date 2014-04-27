package srp.dbgraph;

import java.util.ArrayList;
import java.util.HashMap;

public class DeBruijnGraph {
//	private OpenIntIntHashMap
	private ArrayList<String> allNodes;
	private ArrayList<Integer> allLength;
	private ArrayList<ArrayList<Integer>> allEdges;
	
	private HashMap<Integer, String> allNodesHash;
	private HashMap<Integer, ArrayList<Integer>> allEdgesHash;
//	private int maxNode;
//	private CompatibleSets allSets;
	private int maxNodeIndex;
	private int totalNode;
	
	public DeBruijnGraph() {
		allNodesHash = new HashMap<Integer, String>();
		allEdgesHash = new HashMap<Integer, ArrayList<Integer>>();
		
//		maxNode = 0;
	}
	
	public void addNode(int nodeIndex, String seqs) {
		allNodesHash.put(nodeIndex, seqs);
//		allLength.a(nodeIndex, seqs.length());
		
		allEdgesHash.put(nodeIndex, new ArrayList<Integer>());
		if(nodeIndex> maxNodeIndex){
			maxNodeIndex = nodeIndex;
		}
	}

	public void addEdge(int node, int node2) {
		ArrayList<Integer> arrayList = allEdgesHash.get(node);
		arrayList.add(node2);
		
	}
	
	public void preprocess() {
		checkErrors();
		totalNode = maxNodeIndex+1;
		allNodes = new ArrayList<String>(totalNode);
		allLength = new ArrayList<Integer>(totalNode);
		allEdges = new ArrayList<ArrayList<Integer>>(totalNode);
		
		for (int i = 0; i < totalNode; i++) {
			String string = allNodesHash.get(i);
			allNodes.add(string);
			allLength.add(string.length());
			ArrayList<Integer> arrayList = allEdgesHash.get(i);
			allEdges.add(arrayList);
		}
		
		
	}

	private void checkErrors() {
		if(allNodesHash.size() != (maxNodeIndex+1)){
			throw new IllegalArgumentException("allNodesHash.size() != (maxNodeIndex+1). Skip/miss a node somewhere\t"+allNodes.size() +"\tMax+1:"+ (maxNodeIndex+1) );
		}
//		if(allNodesHash.size() != allLength.size()){
//			throw new IllegalArgumentException("allNodes.size() != allLength.size()\t"+allNodes.size() +"\t"+  allLength.size());
//		}
		if(allNodesHash.size() != allEdgesHash.size()){
			throw new IllegalArgumentException("allNodesHash.size() != allEdgesHash.size()\t"+allNodes.size() +"\t"+ allEdges.size());
		}
		for (ArrayList<Integer> edges : allEdgesHash.values()) {
			for (Integer integer : edges) {
				if(integer > maxNodeIndex){
					throw new IllegalArgumentException("value>maxNodeIndex\t"+integer +"\tMax:"+  maxNodeIndex);
				}
			}
			
		}
	}

	public ArrayList<String> getAllNodes() {
		return allNodes;
	}

	public ArrayList<Integer> getAllLength() {
		return allLength;
	}

	public ArrayList<ArrayList<Integer>> getAllEdges() {
		return allEdges;
	}

	public int getSize() {
		
		return allNodes.size();
	}

	
	
	
}
