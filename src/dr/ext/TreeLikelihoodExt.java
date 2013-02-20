package dr.ext;

import java.util.Arrays;
import java.util.logging.Logger;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AbstractTreeLikelihood;
import dr.evomodel.treelikelihood.AminoAcidLikelihoodCore;
import dr.evomodel.treelikelihood.GeneralLikelihoodCore;
import dr.evomodel.treelikelihood.LikelihoodCore;
import dr.evomodel.treelikelihood.NativeAminoAcidLikelihoodCore;
import dr.evomodel.treelikelihood.NativeGeneralLikelihoodCore;
import dr.evomodel.treelikelihood.NativeNucleotideLikelihoodCore;
import dr.evomodel.treelikelihood.NucleotideLikelihoodCore;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodel.treelikelihood.TreeLikelihood.SiteLikelihoodsStatistic;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;

public class TreeLikelihoodExt extends TreeLikelihood {

	public TreeLikelihoodExt(PatternList patternList, TreeModel treeModel,
			SiteModel siteModel, BranchRateModel branchRateModel,
			TipStatesModel tipStatesModel, boolean useAmbiguities,
			boolean allowMissingTaxa, boolean storePartials,
			boolean forceJavaCore, boolean forceRescaling) {
		super(patternList, treeModel, siteModel, branchRateModel,
				tipStatesModel, useAmbiguities, allowMissingTaxa,
				storePartials, forceJavaCore, forceRescaling);
	}

	/**
	 * 
	 */
	private static final long serialVersionUID = 6348244745369755426L;

	public void updatePatternList(PatternList patternList) {
		

        this.patternList = patternList;
//        this.dataType = patternList.getDataType();
        patternCount = patternList.getPatternCount();
//        stateCount = dataType.getStateCount();

        patternWeights = patternList.getPatternWeights();

//        this.treeModel = treeModel;
        addModel(treeModel);

        nodeCount = treeModel.getNodeCount();

        updateNode = new boolean[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = true;
        }

        likelihoodKnown = false;
		
        patternLogLikelihoods = new double[patternCount];
        resetRootPartials();
//        rootPartials = new double[patternCount * stateCount];

//	}
//
//	public void tTreeLikelihood(PatternList patternList, TreeModel treeModel,
//			SiteModel siteModel, BranchRateModel branchRateModel,
//			TipStatesModel tipStatesModel, 
//			boolean useAmbiguities,
//			boolean allowMissingTaxa, boolean storePartials,
//			boolean forceJavaCore, boolean forceRescaling) {

//		super(TreeLikelihoodParser.TREE_LIKELIHOOD, patternList, treeModel);
		 

//		        super(TreeLikelihoodParser.TREE_LIKELIHOOD, patternList, treeModel);

//		        this.storePartials = storePartials;
		boolean useAmbiguities = false;
		boolean allowMissingTaxa = false;
//		boolean storePartials,
//		boolean forceJavaCore, 
		boolean forceRescaling = false;

		        try {
//		            this.siteModel = siteModel;
//		            addModel(siteModel);

//		            this.frequencyModel = siteModel.getFrequencyModel();
//		            addModel(frequencyModel);

//		            this.tipStatesModel = tipStatesModel;

//		            integrateAcrossCategories = siteModel.integrateAcrossCategories();

		            this.categoryCount = siteModel.getCategoryCount();

//		            final Logger logger = Logger.getLogger("dr.evomodel");
		            String coreName = "Java general";
		            boolean forceJavaCore = false;
		            if (integrateAcrossCategories) {

		                final DataType dataType = patternList.getDataType();

		                if (dataType instanceof dr.evolution.datatype.Nucleotides) {

		                    if (!forceJavaCore && NativeNucleotideLikelihoodCore.isAvailable()) {
		                        coreName = "native nucleotide";
		                        likelihoodCore = new NativeNucleotideLikelihoodCore();
		                    } else {
		                        coreName = "Java nucleotide";
		                        likelihoodCore = new NucleotideLikelihoodCore();
		                    }

		                } 
		            } else {
		                likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
		            }
		//TODO add/remove/change/surpress later            
//		            {
//		              final String id = getId();
//		              logger.info("TreeLikelihood(" + ((id != null) ? id : treeModel.getId()) + ") using " + coreName + " likelihood core");
		//
//		              logger.info("  " + (useAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
//		              logger.info("  With " + patternList.getPatternCount() + " unique site patterns.");
//		            }

//		            if (branchRateModel != null) {
//		                this.branchRateModel = branchRateModel;
////		                logger.info("Branch rate model used: " + branchRateModel.getModelName());
//		            } else {
//		                this.branchRateModel = new DefaultBranchRateModel();
//		            }
//		            addModel(this.branchRateModel);

		            probabilities = new double[stateCount * stateCount];

		            likelihoodCore.initialize(nodeCount, patternCount, categoryCount, integrateAcrossCategories);

		            int extNodeCount = treeModel.getExternalNodeCount();
		            int intNodeCount = treeModel.getInternalNodeCount();
		            /*
		            if (tipStatesModel != null) {
		                tipStatesModel.setTree(treeModel);

		                tipPartials = new double[patternCount * stateCount];

		                for (int i = 0; i < extNodeCount; i++) {
		                    // Find the id of tip i in the patternList
		                    String id = treeModel.getTaxonId(i);
		                    int index = patternList.getTaxonIndex(id);

		                    if (index == -1) {
		                        throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
		                                ", is not found in patternList, " + patternList.getId());
		                    }

		                    tipStatesModel.setStates(patternList, index, i, id);
		                    likelihoodCore.createNodePartials(i);
		                }

		                addModel(tipStatesModel);
		                //useAmbiguities = true;
		            } else */
		            {
		                for (int i = 0; i < extNodeCount; i++) {
		                    // Find the id of tip i in the patternList
		                    String id = treeModel.getTaxonId(i);
		                    int index = patternList.getTaxonIndex(id);

		                    if (index == -1) {
		                        if (!allowMissingTaxa) {
		                            throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
		                                    ", is not found in patternList, " + patternList.getId());
		                        }
		                        if (useAmbiguities) {
		                            setMissingPartials(likelihoodCore, i);
		                        } else {
		                            setMissingStates(likelihoodCore, i);
		                        }
		                    } else {
		                        if (useAmbiguities) {
		                            setPartials(likelihoodCore, patternList, categoryCount, index, i);
		                        } else {
		                            setStates(likelihoodCore, patternList, index, i);
		                        }
		                    }
		                }
		            }
		            for (int i = 0; i < intNodeCount; i++) {
		                likelihoodCore.createNodePartials(extNodeCount + i);
		            }

//		            if (forceRescaling) {
//		                likelihoodCore.setUseScaling(true);
////		                logger.info("  Forcing use of partials rescaling.");
//		            }

		        } catch (TaxonList.MissingTaxonException mte) {
		            throw new RuntimeException(mte.toString());
		        }

		        addStatistic(new SiteLikelihoodsStatistic());
	}



}
