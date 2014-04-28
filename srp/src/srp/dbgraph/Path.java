package srp.dbgraph;

import java.util.ArrayList;

import com.google.common.base.Splitter;

public class Path {
	
	private ArrayList<Integer> nodeList;

	public Path(CompatibleNode[] nodes) {
		nodeList = new ArrayList<Integer>();
	}
	
	public Path(String nodeString){
		nodeList = new ArrayList<Integer>();

//		1:6 7 2:7 30 3:14 42 4:21 88 5:19 182 6:27 241 7:34 258 8:45 317 9:46 328 10:53 374 11:47 448 12:56 545 13:62 612 14:63 701 15:71 792 16:76 798 17:85 831 18:92 957 19:103 1002 20:109 1009 21:118 1068 22:124 1085 23:131 1118 
		Iterable<String> split = Splitter.on(" ").split(nodeString);
		
		for (String string : split) {
			int index = string.indexOf(":"); 
			if(index > 0){

				int node = Integer.parseInt(string.substring(index+1));
				nodeList.add(node);
			}
		}
		
	}
	public ArrayList<Integer> getNodeList(){
		return nodeList;
	}
	public int getNodeCount(){
		return nodeList.size();
	}
}
