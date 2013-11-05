package test.srp.spectrum.likelihood;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.AlignmentUtils;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.likelihood.SpectrumTreeLikelihood;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.inference.model.Parameter;

public class SpectrumTreeLikelihoodTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testSetting() throws Exception {
		// Sub model
		Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25,
				0.25, 0.25 });
		Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100);

		FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
		HKY hky = new HKY(kappa, f);

		// siteModel
		GammaSiteModel siteModel = new GammaSiteModel(hky);
		Parameter mu = new Parameter.Default(
				GammaSiteModelParser.MUTATION_RATE, 1.0, 0,
				Double.POSITIVE_INFINITY);
		siteModel.setMutationRateParameter(mu);

		// treeLikelihood
		Taxon[] taxa = new Taxon[4];
		for (int i = 0; i < taxa.length; i++) {
			taxa[i] = new Taxon("T"+i);
		}
		SimpleAlignment alignment = new SimpleAlignment();
		alignment.addSequence(new Sequence(taxa[0], "R"));
		alignment.addSequence(new Sequence(taxa[1], "K"));
		alignment.addSequence(new Sequence(taxa[2], "C"));
		alignment.addSequence(new Sequence(taxa[3], "C"));
		SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
//		SpectrumAlignmentModel specturmModel = new SpectrumAlignmentModel(shortReads, hapCount)
		TreeModel treeModel = createTreeModel(taxa);
		TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel,
				siteModel, null, null, false, false, true, false, false);
		System.out.println(treeLikelihood.getLogLikelihood());

	}

	@Test
	public void testSpectrumModel() throws Exception {
		// Sub model
		Parameter freqs = new Parameter.Default(new double[] { 0.25, 0.25,
				0.25, 0.25 });
		Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100);

		FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
		HKY hky = new HKY(kappa, f);

		// siteModel
		GammaSiteModel siteModel = new GammaSiteModel(hky);
		Parameter mu = new Parameter.Default(
				GammaSiteModelParser.MUTATION_RATE, 1.0, 0,
				Double.POSITIVE_INFINITY);
		siteModel.setMutationRateParameter(mu);

		// treeLikelihood
		int taxaCount = 4;
		Taxon[] taxa = new Taxon[taxaCount];
		for (int i = 0; i < taxa.length; i++) {
			taxa[i] = new Taxon("taxa_"+i);
		}
		SimpleAlignment alignment = new SimpleAlignment();
		alignment.addSequence(new Sequence(taxa[0], "R"));
		alignment.addSequence(new Sequence(taxa[1], "K"));
		alignment.addSequence(new Sequence(taxa[2], "C"));
		alignment.addSequence(new Sequence(taxa[3], "C"));
		
		
		TreeModel treeModel = createTreeModel(taxa);
		

		String[] seqs = new String[] {
			"A***AAA...",
			".C**AAA...",
			"..G*AAA...",
			"...TTGC..."};
		
		AlignmentMapping aMap = new AlignmentMapping(AlignmentUtils.createAlignment(seqs));
				
		
		SpectrumAlignmentModel specturmModel = new SpectrumAlignmentModel(aMap, taxaCount);
		SpectrumTreeLikelihood treeLikelihood = new SpectrumTreeLikelihood(specturmModel, treeModel,
				siteModel, null, false, false, true, false, false);
		System.out.println(treeLikelihood.getLogLikelihood());

	}
	private static TreeModel createTreeModel(Taxon[] taxa) {

		SimpleNode[] nodes = new SimpleNode[taxa.length*2-1];
		for (int n = 0; n < nodes.length; n++) {
			nodes[n] = new SimpleNode();
		}

		for (int i = 0; i < taxa.length; i++) {
			nodes[i].setTaxon(taxa[i]);
		}
		
		nodes[4].setHeight(0.01);
		nodes[4].addChild(nodes[0]);
		nodes[4].addChild(nodes[1]);

		nodes[5].setHeight(0.01);
		nodes[5].addChild(nodes[2]);
		nodes[5].addChild(nodes[3]);

//		nodes[5].setTaxon(taxa[3]); // gorilla
//
		nodes[6].setHeight(0.02);
		nodes[6].addChild(nodes[4]);
		nodes[6].addChild(nodes[5]);
//
//		nodes[7].setTaxon(taxa[4]); // orangutan
//
//		nodes[8].setHeight(0.069125);
//		nodes[8].addChild(nodes[6]);
//		nodes[8].addChild(nodes[7]);

//		nodes[9].setTaxon(taxa[5]); // siamang

		SimpleNode root = new SimpleNode();
		root = nodes[6];
//		root.setHeight(0.02);
//		root.addChild(nodes[0]);
//		root.addChild(nodes[3]);

		Tree tree = new SimpleTree(root);
		tree.setUnits(Units.Type.YEARS);

		return new TreeModel(tree); // treeModel
	}

}
