package core;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

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
		testReadFiles(dataDir);
		
		
	}

	private static void testReadFiles(String dataDir) {
		
		String trueAlignmentFile = "121101_true_seqs.fasta";
//		String trueAlignmentFile = "zz.fasta";
		String truePhylogenyFile = "121101_true_tree.newick";
		String shortReadFile = "121101_short_reads_10.fasta";
		String refSeqFile = "121101_ref.fasta";
		
		DataImporter dataImporter = new DataImporter(dataDir);
		SimpleAlignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
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
