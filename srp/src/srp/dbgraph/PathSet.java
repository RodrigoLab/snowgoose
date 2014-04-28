package srp.dbgraph;

import java.util.ArrayList;

public class PathSet {
	
	ArrayList<Path> allPaths;
	private boolean offSet;
	
	public PathSet(boolean offSet){
		this.offSet = offSet;
		allPaths  = new ArrayList<Path>();
		if(offSet){
			allPaths.add(null);
		}
	}
	
	public void addPath(Path path){
		allPaths.add(path);
	}
	
	public Path getPath(int i){
		return allPaths.get(i);
	}

	public int getPathCount() {
		return allPaths.size();
	}

	public boolean getOffSet() {
		return offSet;
	}
	
}
