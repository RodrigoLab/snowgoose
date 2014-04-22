package srp.dbgraph;

import java.util.Iterator;

import com.google.common.base.Splitter;

public class CompatibleSet {

	private int node;
	private int depth;

	public CompatibleSet(Iterable<String> split2) {

//		for (String string : split2) {
//			int node = 
//		}
		Iterator<String> iterator = split2.iterator();
		node = Integer.parseInt( iterator.next());
		depth = Integer.parseInt( iterator.next());
		System.out.println(node +"\t"+ depth);
	}

	public void addCompatibleNode(String string) {
		Iterable<String> split2 = Splitter.on(":").split(string);
		Iterator<String> iterator = split2.iterator();
		int cNode = Integer.parseInt( iterator.next());
		int cDepth = Integer.parseInt( iterator.next());
		int cCount = Integer.parseInt( iterator.next());
		System.out.print(cNode +"\t"+ cDepth +"\t"+ cCount +"\t");
		
	}
	
}
