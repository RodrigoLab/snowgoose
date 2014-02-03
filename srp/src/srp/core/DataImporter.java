package srp.core;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import srp.dr.evolution.datatype.ShortReads;
import srp.dr.evolution.io.ShortReadImporter;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.FastaImporter;
import dr.evolution.io.NewickImporter;
//import dr.evolution.io.Importer.ImportException;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.SequenceList;
import dr.evolution.sequence.Sequences;
import dr.evolution.tree.Tree;

public class DataImporter {

//	private SimpleAlignment alignment;
//	private FlexibleTree tree;
	
	
	private String dataDir;
	
	public DataImporter(String dataDir){
		if(dataDir.lastIndexOf(File.separator)!= dataDir.length()){
			dataDir += File.separator;
		}
		this.dataDir = dataDir;
	}
	
//	public DataImporter() {
//		this.dataDir = "";
//	}
	public Alignment importShortReads(String fileName) throws Exception{
		
		Alignment alignment = importShortReads(dataDir, fileName);
		return alignment;
	}
	
	public Alignment importAlignment(String fileName) throws Exception{
		
		Alignment alignment = importAlignment(dataDir, fileName);
		return alignment;
	}
	
	public Tree importTree(String fileName) throws Exception{
		
		Tree tree =  importTree(dataDir, fileName);
		return tree;
	}

	public Sequences importSequence(String fileName) throws Exception{
		
		Sequences seqs =  (Sequences) importSequences(dataDir, fileName);
		return seqs;
	}
	
	public Sequence importRefSeq(String fileName) throws Exception{
		Sequence seq = importRefSeq(dataDir, fileName);
		return seq;
	}
	
	public static Sequence importRefSeq(String dataDir, String fileName) throws Exception{

		
		String filePath = dataDir+fileName;
		
		BufferedReader in = new BufferedReader(new FileReader(filePath));
		FastaImporter importer = new FastaImporter(in, Nucleotides.INSTANCE);
		SequenceList seqList = importer.importSequences();
		Sequence seq = seqList.getSequence(0);
		if(seqList.getSequenceCount()>1){
			System.err.println("Multiple ref seq\t"+seqList.getSequenceCount());
		}
	
		return seq;
	}

	public static SequenceList importSequences(String dataDir, String fileName) throws Exception {

		
		String filePath = dataDir+fileName;
		
		BufferedReader in = new BufferedReader(new FileReader(filePath));
		FastaImporter importer = new FastaImporter(in, Nucleotides.INSTANCE);
		SequenceList seqs = importer.importSequences();

		return seqs;
		
	}

	public static Alignment importAlignment(String dataDir, String fileName) throws Exception{
		
		String filePath = dataDir+fileName;
		
		BufferedReader in = new BufferedReader(new FileReader(filePath));
		FastaImporter importer = new FastaImporter(in, Nucleotides.INSTANCE);
		Alignment alignment = importer.importAlignment();

		return alignment;
	}
	
	public static Alignment importShortReads(String dataDir, String fileName) throws Exception{
		
		String filePath = dataDir+fileName;
		
		BufferedReader in = new BufferedReader(new FileReader(filePath));
		ShortReadImporter importer = new ShortReadImporter(in, ShortReads.INSTANCE);
		Alignment alignment = importer.importAlignment();

		return alignment;
	}
	
		
	public static Tree importTree(String dataDir, String fileName) throws Exception{
		
		String filePath = dataDir+fileName;

		BufferedReader in = new BufferedReader(new FileReader(filePath));
		NewickImporter importer = new NewickImporter(in);
		Tree tree = importer.importTree(null);

		return tree;
	}
	

	
	
}
