package srp.dbgraph;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.StringTokenizer;

import srp.core.DataImporter;

import com.google.common.base.Splitter;

public class DeBruijnImporter extends DataImporter {

	public DeBruijnImporter(String dataDir) {
		super(dataDir);

	}

	
	public DeBruijnGraph importDeBruijnGraph(String fileName) throws Exception{
		DeBruijnGraph dbg = importDeBruijnGraph(dataDir, fileName);
		return dbg;
	}
	
	public CompatibleSets importCompatibleSet(String fileName) throws Exception {
		CompatibleSets allSets = importCompatibleSet(dataDir, fileName);
		return allSets;
		
	}


	public static DeBruijnGraph importDeBruijnGraph(String dataDir, String fileName) throws Exception{
		
		
		DeBruijnGraph dbg = new DeBruijnGraph();
		BufferedReader inFile  = new BufferedReader(new FileReader(dataDir+fileName));
		String inline;
		int length = 0;
//		SpectrumAlignmentModel spectrumModel = null;
		
		String regex = "\\d+\\s\\d+";
		while((inline = inFile.readLine())!=null){
			StringTokenizer st = new StringTokenizer(inline);
			if(inline.matches(regex)){
//				System.out.println("Secnod half\t"+inline);
//				StringTokenizer st = new StringTokenizer(inline);
				int node1 = Integer.parseInt(st.nextToken());
				int node2 = Integer.parseInt(st.nextToken());
				dbg.addEdge(node1, node2);
				
			}
			else{
				int node = Integer.parseInt(st.nextToken());
				String seqs = st.nextToken();
				dbg.addNode(node, seqs);
			}
		}
//		dbg.test();
		inFile.close();
		return dbg;
		
	}

	
	public static CompatibleSets importCompatibleSet(String dataDir, String fileName) throws Exception{
//		DeBruijnGraph dbg = new DeBruijnGraph();
//		CompatibleSet compSet = new CompatibleSet();
//		ArrayList<CompatibleNode> allSets = new ArrayList<CompatibleNode>();
		CompatibleSets compatibleSets = new CompatibleSets();
		BufferedReader inFile  = new BufferedReader(new FileReader(dataDir+fileName));
		String inline;
		int length = 0;
//		SpectrumAlignmentModel spectrumModel = null;
		
		String regex = "\\d+\\s\\d+";
		while((inline = inFile.readLine())!=null){
			
			Iterable<String> split = Splitter.on(" ").split(inline.trim());
			CompatibleNode compNode = null;//new CompatibleSet();
			for (String string : split) {
				if(compNode == null){
					Iterable<String> split2 = Splitter.on(":").split(string);
					compNode = new CompatibleNode(split2);
				}
				else{
					compNode.addCompatibleNode(string);
				}
			}
			compatibleSets.addCompatibleNode(compNode);
			
//			4:30 4:0:4903 18:31:554 30:32:20056 17:61:3745 27:66:28245 41:126:106 55:128:55372 63:405:3356 72:429:218 82:431:3350 

//			11:31 11:0:4060 18:31:462 29:32:97682 
//			<node no>:<depth node> 
//			<compatible node>:<node depth>:<count node> 
//			<compatible node>:<node depth>:<count node> <compatible node>:<node depth>:<count node>  ...

			
		}
		inFile.close();
		
		return compatibleSets;
	}

	public static PathSet importPathSets(String dataDir, String fileName, boolean offSet) throws Exception{

		BufferedReader inFile  = new BufferedReader(new FileReader(dataDir+fileName));
		PathSet pathSet = new PathSet(offSet);
		String inline;
		
		while((inline = inFile.readLine())!=null){
			if(!inline.startsWith(">")){
				Path p = new Path(inline);
				pathSet.addPath(p);
			}
		}
		inFile.close();
		return pathSet;
	}


	public PathSet importPathSets(String fileName, boolean b) throws Exception {
		
		return importPathSets(dataDir, fileName, b);
	}
}
