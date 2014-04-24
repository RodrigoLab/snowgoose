package srp.dbgraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import com.google.common.base.Splitter;

public class CompatibleNode {

	private int nodeIndex;
	private int nodeDepth;
	private ArrayList<Integer> cNodeList;
	private HashMap<Integer, Integer> cCount;
	private HashMap<Integer, Integer> cDepth;
	
	public CompatibleNode(Iterable<String> split2) {


		Iterator<String> iterator = split2.iterator();
		nodeIndex = Integer.parseInt( iterator.next());
		nodeDepth = Integer.parseInt( iterator.next());
		cNodeList = new ArrayList<>();
		cCount = new HashMap<Integer, Integer>();
		cDepth= new HashMap<Integer, Integer>();
		System.out.println(nodeIndex +"\t"+ nodeDepth);
	}

	public void addCompatibleNode(String string) {
		Iterable<String> split2 = Splitter.on(":").split(string);
		Iterator<String> iterator = split2.iterator();
		int cNode = Integer.parseInt( iterator.next());
		int tempDepth = Integer.parseInt( iterator.next());
		int tempCount = Integer.parseInt( iterator.next());
//		System.out.print(cNode +"\t"+ cDepth +"\t"+ cCount +"\t");
//		int[] cDCArray = new int[]{cDepth, cCount};
//		compatible.put(cNode, cDCArray);
		cNodeList.add(cNode);
		cDepth.put(cNode, tempDepth);
		cCount.put(cNode, tempCount);
		
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

	public ArrayList<Integer> getCNodeList(){
		return cNodeList;
	}
	public int getNodeIndex() {
		return nodeIndex;
	}

	public int getDepth() {
		return nodeDepth;
	}
	
}
