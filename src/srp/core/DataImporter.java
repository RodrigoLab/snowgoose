package srp.core;

import java.io.BufferedReader;
import java.io.EOFException;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.FastaImporter;
import dr.evolution.io.NewickImporter;
//import dr.evolution.io.Importer.ImportException;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.SequenceList;
import dr.evolution.sequence.Sequences;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;

public class DataImporter {

//	private SimpleAlignment alignment;
//	private FlexibleTree tree;
	
	
	private String dataDir;
	
	public DataImporter(String dataDir){
		this.dataDir = dataDir;
	}
	
	public DataImporter() {
		this.dataDir = "";
	}

	public Alignment importAlignment(String fileName){
		
		Alignment alignment = importAlignment(dataDir, fileName);
		return alignment;
	}
	
	public Tree importTree(String fileName){
		
		Tree tree =  importTree(dataDir, fileName);
		return tree;
	}

	public Sequences importSequence(String fileName){
		
		Sequences seqs =  (Sequences) importSequences(dataDir, fileName);
		return seqs;
	}
	
	public Sequence importRefSeq(String fileName){
		Sequence seq = importRefSeq(dataDir, fileName);
		return seq;
	}
	
	public static Sequence importRefSeq(String dataDir, String fileName) {

		Sequence seq = null;
		String filePath = dataDir+fileName;
		try {
			BufferedReader in = new BufferedReader(new FileReader(filePath));
			FastaImporter importer = new FastaImporter(in, Nucleotides.INSTANCE);
			SequenceList seqList = importer.importSequences();
			seq = seqList.getSequence(0);
			if(seqList.getSequenceCount()>1){
				System.err.println("Multiple ref seq\t"+seqList.getSequenceCount());
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return seq;
	}

	public static SequenceList importSequences(String dataDir, String fileName) {

		SequenceList seqs = null;
		String filePath = dataDir+fileName;
		try {
			BufferedReader in = new BufferedReader(new FileReader(filePath));
			FastaImporter importer = new FastaImporter(in, Nucleotides.INSTANCE);
			seqs = importer.importSequences();

		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return seqs;
		
	}

	public static Alignment importAlignment(String dataDir, String fileName){
		
		Alignment alignment = null;
		String filePath = dataDir+fileName;
		try {
			BufferedReader in = new BufferedReader(new FileReader(filePath));
			FastaImporter importer = new FastaImporter(in, Nucleotides.INSTANCE);
			alignment =  importer.importAlignment();

		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return alignment;
	}
	
		
	public static Tree importTree(String dataDir, String fileName){
		
		Tree tree = null;
		String filePath = dataDir+fileName;
		try {
			BufferedReader in = new BufferedReader(new FileReader(filePath));
			NewickImporter importer = new NewickImporter(in);
			tree = importer.importTree(null);

//			for (int i = 0; i < tree.getNodeCount(); i++) {
//				System.out.println("Node\t"+i+"\t"+tree.getBranchLength(tree.getNode(i)));
//			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return tree;
	}
	

//			List<Taxon> taxonList = alignment.asList();			
//			for (int i = 0; i < alignment.getSequenceCount(); i++) {
//				System.out.println(alignment.getTaxon(i).toString());
//				System.out.println(alignment.getSequence(i).getSequenceString());
//				
//			}
//			for (Taxon taxon : taxonList) {
//				System.out.println(alignment.getTaxonIndex(taxon));
//			}
//			
//			List<Taxon> treeTaxon = tree.asList();
//			for (Taxon taxon : treeTaxon) {
//				System.out.println(taxon.toString());
//				
//			}
//			System.out.println(tree.getInternalNodeCount());
//			System.out.println(tree.getExternalNodeCount());
//			System.out.println(tree.getRoot());
//			System.out.println(tree.toString());
			


//	/**
//	 * @return the tree
//	 */
//	public FlexibleTree getTree() {
//		return tree;
//	}
//
//	public SimpleAlignment getAlignment() {
//		return alignment;
//	}
//	
	
	
	
}
