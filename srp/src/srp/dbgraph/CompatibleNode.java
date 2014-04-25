package srp.dbgraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import cern.colt.map.OpenIntIntHashMap;

import com.google.common.base.Splitter;
import com.google.common.collect.Collections2;
import com.google.common.primitives.Ints;
import com.google.common.primitives.Primitives;

public class CompatibleNode {

	private int nodeIndex;
	private int nodeDepth;
	private ArrayList<Integer> cNodeList;
//	private HashMap<Integer, Integer> cCount;
//	private HashMap<Integer, Integer> cDepth;
	
//	private OpenIntIntHashMap cCountMap;
	private OpenIntIntHashMap cCount;
	private OpenIntIntHashMap cDepth;
	
	
	private int[] cCountArray;
	private int[] cDepthArray;
	private int[] cNodeArray;
	
			
	
	public CompatibleNode(Iterable<String> split2) {


		Iterator<String> iterator = split2.iterator();
		nodeIndex = Integer.parseInt( iterator.next());
		nodeDepth = Integer.parseInt( iterator.next());
		cNodeList = new ArrayList<>();
//		cCount = new HashMap<Integer, Integer>();
//		cDepth= new HashMap<Integer, Integer>();
		cCount = new OpenIntIntHashMap();
		cDepth = new OpenIntIntHashMap();
		
		System.out.println(nodeIndex +"\t"+ nodeDepth);
	}

	public void addCompatibleNode(String string) {
		Iterable<String> split2 = Splitter.on(":").split(string);
		Iterator<String> iterator = split2.iterator();
		int cNode = Integer.parseInt( iterator.next());
		int tempDepth = Integer.parseInt( iterator.next());
		int tempCount = Integer.parseInt( iterator.next());

		cNodeList.add(cNode);
		cDepth.put(cNode, tempDepth);
		cCount.put(cNode, tempCount);
		
//		cCountMap.put(cNode, tempCount);
	}
	
	public void proprocess() {
		int totalCNodeCount = cNodeList.size();
		if( cNodeList.size() != cDepth.size() ){
			throw new IllegalArgumentException("cNodeList.size() != cDepth.size()\t"+ cNodeList.size() +"\t"+ cDepth.size());
		}
		if( cDepth.size() != cCount.size() ){
			throw new IllegalArgumentException("cDepth.size() != cCount.size()\t"+ cDepth.size() +"\t"+  cCount.size());
		}
		
//		cNodeArray = new int[totalCNodeCount];
//		cCountArray = new int[totalCNodeCount];
//		cDepthArray = new int[totalCNodeCount];
//		cNodeArray = cNodeList.toArray(cNodeArray);
		cNodeArray = Ints.toArray(cNodeList);
	}
	
	public int getCNodeDepth(int cNode){
		return cDepth.get(cNode);
	}
	public int getCNodeCount(int cNode){
		return cCount.get(cNode);
	}
//	public int[] getCNodeInfo(int cNode){
//		return compatible.get(cNode);
//	}

	@Deprecated
	public ArrayList<Integer> getCNodeList(){
		return cNodeList;
	}
	
	public int[] getCNodeArray(){
		return cNodeArray;
	}
	
	public int getNodeIndex() {
		return nodeIndex;
	}

	public int getDepth() {
		return nodeDepth;
	}
	
}
