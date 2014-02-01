package test.dr.rj;


import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import com.google.common.base.Strings;

import srp.core.DataImporter;
import srp.haplotypes.HaplotypeModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.rj.SSTreeModel;

public class SSTreeModelTest {
	
	private static Tree truePhylogeny;
	private static SimpleTree tree;
	private static Alignment shortReads;
	private SSTreeModel treeModel;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		

		String dataDir = "/home/sw167/workspaceSrp/ABI/unittest/";
		int noOfRecoveredHaplotype=6;
		String shortReadFile = "H6_srp.fasta";
		String treeFile = "H6_haplotypes.tree";
				
		DataImporter dataImporter = new DataImporter(dataDir);
		truePhylogeny = dataImporter.importTree(treeFile);
		TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, truePhylogeny, false);
		
		shortReads = dataImporter.importAlignment(shortReadFile);
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
//		treeModel = new RJTreeModel(tree);// treeModel
		treeModel = new SSTreeModel(TreeModel.TREE_MODEL, truePhylogeny, false);
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testConstructorFile() throws Exception {
		treeModel = new SSTreeModel(TreeModel.TREE_MODEL, truePhylogeny, false);
		NodeRef root = treeModel.getRoot();
		NodeRef trueRoot = treeModel.getTrueRoot();
		if(root.getNumber() != trueRoot.getNumber()){
			assertEquals(treeModel.getChild(root, 0).getNumber(), trueRoot.getNumber());
		}
        for (int t = 0; t < treeModel.getInternalNodeCount(); t++) {

            NodeRef node =  treeModel.getInternalNode(t);
            assertEquals(treeModel.getChildCount(node), 2);

        }
        assertEquals(6, treeModel.getTrueExternalNodeCount());
        assertEquals(5, treeModel.getTrueInternalNodeCount());
        assertEquals(16, treeModel.getExternalNodeCount());
        assertEquals(15, treeModel.getInternalNodeCount());
        
		
//		treeModel.insertNode();

	}
	
	@Test
	public void testConstructor() throws Exception {
		
		for (int i = 0; i < 100; i++) {
			int noOfRecoveredHaplotype;
			do{
				noOfRecoveredHaplotype = MathUtils.nextInt(150);
			}while(noOfRecoveredHaplotype<3);
			tree = createTree(noOfRecoveredHaplotype);
			treeModel = new SSTreeModel(tree);// treeModel

			NodeRef root = treeModel.getRoot();
			NodeRef trueRoot = treeModel.getTrueRoot();
			if(root.getNumber() != trueRoot.getNumber()){
				assertEquals(treeModel.getChild(root, 0).getNumber(), trueRoot.getNumber());
			}
	        for (int t = 0; t < treeModel.getInternalNodeCount(); t++) {
	
	            NodeRef node =  treeModel.getInternalNode(t);
	            assertEquals(treeModel.getChildCount(node), 2);
	
	        }
	        for (int t = 0; t < treeModel.getExternalNodeCount(); t++) {
	        	if(t<noOfRecoveredHaplotype){
	        		assertTrue(treeModel.getTaxonId(t).startsWith("hap_"));
	        	}
	        	else{
	        		assertEquals(Integer.toString(t), treeModel.getTaxonId(t));
	        	}
	        }
	        assertEquals(noOfRecoveredHaplotype, treeModel.getTrueExternalNodeCount());
	        assertEquals(noOfRecoveredHaplotype-1, treeModel.getTrueInternalNodeCount());
   
		}
	}
	
	@Test
	public void testConstructorNoRescale() throws Exception {

		for (int i = 0; i < 100; i++) {
			int power;
			do{
				power = MathUtils.nextInt(8);
			}while(power<2);
			int noOfRecoveredHaplotype = (int) (Math.pow(2, power));
			tree = createTree(noOfRecoveredHaplotype);
			treeModel = new SSTreeModel(tree);

			NodeRef root = treeModel.getRoot();
			NodeRef trueRoot = treeModel.getTrueRoot();

	        for (int t = 0; t < treeModel.getInternalNodeCount(); t++) {
	            NodeRef node =  treeModel.getInternalNode(t);
	            assertEquals(treeModel.getChildCount(node), 2);
	        }
			if(noOfRecoveredHaplotype>15){
				assertEquals(root.getNumber(), trueRoot.getNumber());
		        assertEquals(noOfRecoveredHaplotype, treeModel.getExternalNodeCount());
		        assertEquals(noOfRecoveredHaplotype-1, treeModel.getInternalNodeCount());
		        assertEquals(noOfRecoveredHaplotype, treeModel.getTrueExternalNodeCount());
		        assertEquals(noOfRecoveredHaplotype-1, treeModel.getTrueInternalNodeCount());
			}
			else{
				assertEquals(30, root.getNumber());
				assertEquals(15+noOfRecoveredHaplotype-1, trueRoot.getNumber());
		        assertEquals(16, treeModel.getExternalNodeCount());
		        assertEquals(15, treeModel.getInternalNodeCount());
		        assertEquals(noOfRecoveredHaplotype, treeModel.getTrueExternalNodeCount());
		        assertEquals(noOfRecoveredHaplotype-1, treeModel.getTrueInternalNodeCount());

			}
//	        assertEquals(noOfRecoveredHaplotype, treeModel.getExternalNodeCount());
//	        assertEquals(noOfRecoveredHaplotype-1, treeModel.getInternalNodeCount());
//	        assertEquals(noOfRecoveredHaplotype, treeModel.getTrueExternalNodeCount());
//	        assertEquals(noOfRecoveredHaplotype-1, treeModel.getTrueInternalNodeCount());
	        
		}
	}
	
	@Test
	public void testInsertNode() throws Exception {
		System.out.println(treeModel.getExternalNodeCount() +"\t"+ treeModel.getTrueExternalNodeCount() +"\t"+ treeModel.getNodeCount());
		System.out.println(treeModel.toString());
//		treeModel = treeModel.insertNode();
		treeModel.insertNode();
//		treeModel = treeModel.insertNode();
////		treeModel.insertNode();
//		treeModel.insertNode();
//		treeModel.insertNode();
		System.out.println(treeModel.getExternalNodeCount() +"\t"+ treeModel.getTrueExternalNodeCount() +"\t"+ treeModel.getNodeCount());
		System.out.println(treeModel.toString());
	}

	private static SimpleTree createTree(int noOfRecoveredHaplotype) {
		HaplotypeModel haplotypeModel = new HaplotypeModel(shortReads, noOfRecoveredHaplotype);
		Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000.0, 100, 100000.0);
		ConstantPopulationModel popModel = new ConstantPopulationModel(popSize, Units.Type.YEARS);
		ConstantPopulation constant = (ConstantPopulation) popModel.getDemographicFunction();
		CoalescentSimulator simulator = new CoalescentSimulator();
		SimpleTree simulateTree = simulator.simulateTree(haplotypeModel, constant);
		return simulateTree;
	}
}
