package jars;

import java.io.OutputStream;
import java.io.PrintStream;

import likelihood.LikelihoodCalculation;
import core.DataImporter;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.sequence.Sequences;
import dr.evomodel.tree.TreeModel;

public class MainShortReadLikelihood {

	public static void main(String[] args) {

		if (args.length == 0) {
			System.out
					.println("Usage: java -jar testShortReadLikelihood.jar alignmentFile shortReadFile");
		}

		System.setErr(new PrintStream(new OutputStream() {
			@Override
			public void write(int b) {
			}
		}));

		String trueAlignmentFile = args[0];
		String shortReadFile = args[1];

		DataImporter dataImporter = new DataImporter();
		Alignment trueAlignment = dataImporter.importAlignment(trueAlignmentFile);
		Sequences shortReads = dataImporter.importSequence(shortReadFile);

		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL);

		LikelihoodCalculation li = new LikelihoodCalculation(treeModel, trueAlignment, shortReads);

		System.out.println(li.getShortReadLikelihood());

	}

}
