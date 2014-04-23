package srp.dbgraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import com.google.common.base.Splitter;

public class CompatibleSet {

	private int node;
	private int depth;
	private ArrayList<Integer> cNodeList;
	private HashMap<Integer, Integer> cCount;
	private HashMap<Integer, Integer> cDepth;
	
	public CompatibleSet(Iterable<String> split2) {

//		for (String string : split2) {
//			int node = 
//		}
		Iterator<String> iterator = split2.iterator();
		node = Integer.parseInt( iterator.next());
		depth = Integer.parseInt( iterator.next());
		cNodeList = new ArrayList<>();
		cCount = new HashMap<Integer, Integer>();
		cDepth= new HashMap<Integer, Integer>();
		System.out.println(node +"\t"+ depth);
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
	public int getNode() {
		return node;
	}

	public int getDepth() {
		return depth;
	}
	
}
