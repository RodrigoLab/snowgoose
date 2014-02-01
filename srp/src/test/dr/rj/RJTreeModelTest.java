package test.dr.rj;


import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.core.DataImporter;
import srp.haplotypes.HaplotypeModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.rj.RJTreeModel;
import dr.rj.RJTreeModel.Node;;

public class RJTreeModelTest {

	private static Tree tree;
	private RJTreeModel treeModel;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {

		String dataDir = "/home/sw167/workspaceSrp/ABI/unittest/";
		int noOfRecoveredHaplotype=6;
		String shortReadFile = "H6_srp.fasta";

		DataImporter dataImporter = new DataImporter(dataDir);
		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		
		HaplotypeModel haplotypeModel = new HaplotypeModel(shortReads, noOfRecoveredHaplotype);
		
		// coalescent
		Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000.0, 100, 100000.0);

		// Random treeModel
		ConstantPopulationModel popModel = new ConstantPopulationModel(popSize, Units.Type.YEARS);
		ConstantPopulation constant = (ConstantPopulation) popModel.getDemographicFunction();
		CoalescentSimulator simulator = new CoalescentSimulator();
		tree = simulator.simulateTree(haplotypeModel, constant);
		
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		treeModel = new RJTreeModel(tree);// treeModel
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void testAddTaxon() throws Exception {
		
		System.out.println(treeModel.getExternalNodeCount() +"\t"+ treeModel.getNodeCount());
		treeModel = treeModel.insertNode();
		System.out.println(treeModel.getExternalNodeCount() +"\t"+ treeModel.getNodeCount());
	}

}
