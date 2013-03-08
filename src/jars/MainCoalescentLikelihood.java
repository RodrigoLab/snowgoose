package jars;

import srp.core.DataImporter;
import srp.likelihood.LikelihoodCalculation;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.inference.model.Parameter;
public class MainCoalescentLikelihood {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		if (args.length == 0) {
			System.out
			.println("Usage: java -jar MainCoalescentLikelihood.jar treeFile popSize");
		}

		String truePhylogenyFile = args[0];
		int p = -1;
		try {
			p = Integer.parseInt(args[1]);
		} catch (Exception e) {
			e.printStackTrace();
		}

		DataImporter dataImporter = new DataImporter();
		Tree truePhylogeny = dataImporter.importTree(truePhylogenyFile);

		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false, false);
		Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, p, 0, p*100);

		LikelihoodCalculation li = new LikelihoodCalculation(treeModel, popSize);

		System.out.println(li.getCoalescentLikelhood());
		

	}

}
