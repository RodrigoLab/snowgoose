package core;
/*
 -XX:CompileThreshold=50000 -XX:+CITime
 
*/
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;
import java.util.LinkedHashMap;

import alignment.AlignmentMapping;
import alignment.AlignmentMatrix;

import likelihood.LikelihoodCalculation;

import io.ShortReadImporter;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.FastaImporter;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.Sequences;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;


public class Main {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/Stage0/";
//		System.out.println(System.getProperty("java.library.path"));
//		System.setProperty("java.library.path","/home/sw167/PostdocLarge/Software/BEAST/BEASTv1.7.1/lib/");
//		System.out.println(System.getProperty("java.library.path"));
//		
//		testReadFiles(dataDir);
		
		testAlignmentMapping();
		
//		testAlignmentMappingLikelihood();
		
//		
//		
//		String s = dataDir;
//		char[] c = dataDir.toCharArray();
//		
//		long time1, time2;
//		time1 = System.currentTimeMillis();
//		for (int i = 0; i < 1e8; i++) {
////			for (int j = 0; j < c.length; j++) {
////				int t = c[j]+1;
////			}\
//			for (int j = 0; j < s.length(); j++) {
//				char t = s.charAt(j);
//			}
////			int length = c.length;
//		}
//		time2 = System.currentTimeMillis();
//		System.out.println("\t" + (time2 - time1) + "\t");

	
		
	}
	private static void testAlignmentMappingLikelihood(){
		
		String dataDir = "/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/";

		String trueAlignmentFile = "1110_10_org.phyml";
//		String trueAlignmentFile = "zz.fasta";
		String phylogenyFile = "1110_10_org.phyml_phyml_tree.txt";
		String shortReadFile = "1110_10_align_2.fasta";
		String refSeqFile = "1110_10.ref";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		Tree truePhylogeny = dataImporter.importTree(phylogenyFile);
		Sequences shortReads = dataImporter.importSequence(shortReadFile);
		Sequence refSeq = dataImporter.importRefSeq(refSeqFile);
		
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);
		
		
		LikelihoodCalculation li = new LikelihoodCalculation(treeModel, trueAlignment, shortReads);
		li.setPopSize(3000,0,30000);
//		li.setTreeAndAlignment(treeModel, trueAlignment);
		System.out.println(li.getTreeLikelihood());
		System.out.println(li.getCoalescentLikelhood());

		System.out.println(li.getShortReadLikelihood());
		
		
		
	}

	private static void testAlignmentMapping() {
		String dir = "/home/sw167/Postdoc/Project_A2BI_temp/data/srAlignment/";
		String refFileName = "1110_10.ref";
		String srpFileName = "1110_10_align_test.fasta";
		
		try {
			

			DataImporter dataImporter = new DataImporter(dir);
//			Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
//			Tree truePhylogeny = dataImporter.importTree(phylogenyFile);
			Alignment shortReads = dataImporter.importAlignment(srpFileName);
			Sequence refSeq = dataImporter.importRefSeq(refFileName);
			System.out.println("===");
			System.out.println(shortReads.getPatternCount());
			System.out.println(shortReads.getSiteCount());
			System.out.println(shortReads.getStateCount());
			System.out.println(shortReads.getTaxonCount());
//			System.out.println(shortReads.get
			for (int i = 0; i < shortReads.getSequenceCount(); i++) {
				Sequence s = shortReads.getSequence(i);
				Taxon t = s.getTaxon();
//				System.out.println(s.getId() +"\t"+ t.getId());	
//				System.out.println(s.getChar(0));
//				System.out.println(s.getLength());
				
				
				
//				System.out.println(s.getLength() +"\t"+ s.getSequenceString());
			}
//			BufferedReader in = new BufferedReader(new FileReader(dir+refFileName));
//			in.readLine();
////			String refSeq = in.readLine();
//			
//			in = new BufferedReader(new FileReader(dir+srpFileName));
//			String input = in.readLine();
//			in.readLine();
//			int count = 0;
//			
//			LinkedHashMap<String, String> alignment = new LinkedHashMap();
//			while((input=in.readLine())!= null){
//				if(input.indexOf(">")!= -1){
//					alignment.put(input.substring(1), in.readLine());
//				}
//				count++;
//				if(count==9){
//					break;
//				}
//			}
//			
//			System.out.println(alignment.toString() );
			
			AlignmentMapping amap = new AlignmentMapping(shortReads);
//			AlignmentMapping amap = new AlignmentMapping(1200);
//			for (int i = 0; i < shortReads.getSequenceCount(); i++) {
//				Sequence s = shortReads.getSequence(i);
//				amap.addSequence(s);
//			}
//			shortReads.
//			
//			for (String key: alignment.keySet()){
////				System.out.println(key);
//				
////				System.out.println(index);// +"\t"+ alignment.get(key));
//				amap.addSeq(key, alignment.get(key));
//			}
//			System.out.println(am.toString());
			AlignmentMatrix alignmentMatrix = new AlignmentMatrix(amap, 5);
			alignmentMatrix.testGetSeq();
			
			alignmentMatrix.testMultipleSeq();
			alignmentMatrix.toAlignment();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}

	private static void testReadFiles(String dataDir) {
		
		String trueAlignmentFile = "121101_true_seqs.fasta";
//		String trueAlignmentFile = "zz.fasta";
		String truePhylogenyFile = "121101_true_tree.newick";
		String shortReadFile = "121101_short_reads_10.fasta";
		String refSeqFile = "121101_ref.fasta";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		Tree truePhylogeny = dataImporter.importTree(truePhylogenyFile);
		Sequences shortReads = dataImporter.importSequence(shortReadFile);
		Sequence refSeq = dataImporter.importRefSeq(refSeqFile);
		
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);
		
		
		LikelihoodCalculation li = new LikelihoodCalculation(treeModel, trueAlignment, shortReads);
		li.setPopSize(3000,0,30000);
//		li.setTreeAndAlignment(treeModel, trueAlignment);
		System.out.println(li.getTreeLikelihood());
		System.out.println(li.getCoalescentLikelhood());

		System.out.println(li.getShortReadLikelihood());
		
	}

}
