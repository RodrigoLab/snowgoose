/*
 * TreeLikelihood.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package srp.likelihood.spectrum.treelikelihood;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperationRecord;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evolution.util.TaxonList.MissingTaxonException;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.LikelihoodCore;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inference.model.Model;

/**
 * SpectrumTreeLikelihood - implements a Likelihood Function for specturms on a tree.
 *
 * @author Steven Wu
 */

public class SpectrumTreeLikelihood extends AbstractSpectrumTreeLikelihood {
    /**
	 * 
	 */
	private static final long serialVersionUID = -5701472737013978020L;

	private static final boolean DEBUG = false;

	private int updateExternalNodeIndex;

    /**
     * Constructor.
     */
    public SpectrumTreeLikelihood(
//    		PatternList patternList,
    		SpectrumAlignmentModel spectrumModel,
                          TreeModel treeModel,
                          SiteModel siteModel,
                          BranchRateModel branchRateModel,
//                          TipStatesModel tipStatesModel,
                          boolean useAmbiguities,
                          boolean allowMissingTaxa,
                          boolean storePartials,
                          boolean forceJavaCore,
                          boolean forceRescaling) {

        super(TreeLikelihoodParser.TREE_LIKELIHOOD, spectrumModel, treeModel);

        this.storePartials = storePartials;
        updateExternalNodeIndex = -1;
        try {
        	
            this.siteModel = siteModel;
            addModel(siteModel);

            this.frequencyModel = siteModel.getFrequencyModel();
            addModel(frequencyModel);

            this.categoryCount = siteModel.getCategoryCount();

            final Logger logger = Logger.getLogger("dr.evomodel");
            logger.setLevel(Level.OFF);
            final DataType dataType = spectrumModel.getDataType();
            
            String coreName = "Java general";
            coreName = "Java spectrumnucleotide";
            likelihoodCore = new SpectrumNucleotideLikelihoodCore();
//            if (integrateAcrossCategories) {
//
//                final DataType dataType = patternList.getDataType();
//
//                if (dataType instanceof dr.evolution.datatype.Nucleotides) {
//
//                    if (!forceJavaCore && NativeNucleotideLikelihoodCore.isAvailable()) {
//                        coreName = "native nucleotide";
//                        likelihoodCore = new NativeNucleotideLikelihoodCore();
//                    } else {
//                        coreName = "Java nucleotide";
//                        likelihoodCore = new NucleotideLikelihoodCore();
//                    }
//
            {
              final String id = getId();
              logger.info("TreeLikelihood(" + ((id != null) ? id : treeModel.getId()) + ") using " + coreName + " likelihood core");

              logger.info("  " + (useAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
//              logger.info("  With " + spectrumModel.getPatternCount() + " unique site patterns.");
              logger.info("  With " + spectrumModel.getSpectrumLength() + " sites.");
            }

            if (branchRateModel != null) {
                this.branchRateModel = branchRateModel;
                logger.info("Branch rate model used: " + branchRateModel.getModelName());
            } else {
                this.branchRateModel = new DefaultBranchRateModel();
            }
            addModel(this.branchRateModel);

            probabilities = new double[stateCount * stateCount];

            likelihoodCore.initialize(nodeCount, patternCount, categoryCount, integrateAcrossCategories);

            int extNodeCount = treeModel.getExternalNodeCount();
            int intNodeCount = treeModel.getInternalNodeCount();

            for (int i = 0; i < extNodeCount; i++) {
                // Find the id of tip i in the patternList
                String id = treeModel.getTaxonId(i);
                int index = spectrumModel.getTaxonIndex(id);
                if (index == -1) {
                      throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
                              ", is not found in patternList, " + spectrumModel.getId());
                }
                setPartials(likelihoodCore, spectrumModel, index, i);
            }

            for (int i = 0; i < intNodeCount; i++) {
                likelihoodCore.createNodePartials(extNodeCount + i);
            }

            if (forceRescaling) {
                likelihoodCore.setUseScaling(true);
                logger.info("  Forcing use of partials rescaling.");
            }
            
        } catch (TaxonList.MissingTaxonException mte) {
            throw new RuntimeException(mte.toString());
        }

        
//        addStatistic(new SiteLikelihoodsStatistic());
    }

    public final LikelihoodCore getLikelihoodCore() {
        return likelihoodCore;
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    /**
     * Handles model changed events from the submodels.
     */
    @Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
//    	System.out.println(model.getModelName());
//    	System.out.println(object.getClass());
        if (model == treeModel) {
            if (object instanceof TreeModel.TreeChangedEvent) {

                if (((TreeModel.TreeChangedEvent) object).isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    updateNodeAndChildren(((TreeModel.TreeChangedEvent) object).getNode());

                } else if (((TreeModel.TreeChangedEvent) object).isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    updateAllNodes();
                } else {
                    // Other event types are ignored (probably trait changes).
                    //System.err.println("Another tree event has occured (possibly a trait change).");
                }
            }

        } else if (model == branchRateModel) {
            if (index == -1) {
                updateAllNodes();
            } else {
                if (DEBUG) {
                if (index >= treeModel.getNodeCount()) {
                    throw new IllegalArgumentException("Node index out of bounds");
                }
                }
                updateNode(treeModel.getNode(index));
            }

        } else if (model == frequencyModel) {

            updateAllNodes();


        } else if (model instanceof SiteModel) {

            updateAllNodes();

        } 

        else if (model == spectrumModel){
        	SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
			int spectrumIndex = record.getSpectrumIndex();
//			int siteIndex = record.getAllSiteIndexs()[0];
//        	updateExternalNodeIndex = -1;
//        	spectrumIndex -> taxonName -> indexOnTree
			String taxonId = spectrumModel.getTaxonId(spectrumIndex);
            updateExternalNodeIndex = treeModel.getTaxonIndex(taxonId );
//            int index = spectrumModel.getTaxonIndex(id);

            if (updateExternalNodeIndex == -1) {
            	try {
					throw new TaxonList.MissingTaxonException("Taxon, " + taxonId + ", in tree, " + treeModel.getId() +
					          ", is not found in patternList, " + spectrumModel.getId());
				} catch (MissingTaxonException e) {
					e.printStackTrace();
				}
            }
//            updateAllNodes();
            double[] partials = new double[patternCount * stateCount];
            int v = 0;
            //TODO only update one part of the array
//            System.out.println("Update Node: "+updateExternalNodeIndex +"\tSpectrumIndex: "+ spectrumIndex);
            for (int i = 0; i < patternCount; i++) {
                double[] frequencies = spectrumModel.getSpecturmFrequencies(spectrumIndex, i);
                //TODO use siteIndex here
//                System.out.println(Arrays.toString(frequencies));
                for (int j = 0; j < stateCount; j++) {
                	partials[v] = frequencies[j];
                    v++;
                }
            }
//            updateAllNodes();
//            updateAllPatterns();
//            makeDirty();
            likelihoodCore.setNodePartialsForUpdate(updateExternalNodeIndex);
            likelihoodCore.setCurrentNodePartials(updateExternalNodeIndex, partials);

            
        }
        
        
        
        
        else {

            throw new RuntimeException("Unknown componentChangedEvent");
        }

        super.handleModelChangedEvent(model, object, index);
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the additional state other than model components
     */
    @Override
	protected void acceptState() {
    	SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();

	}
    @Override
	protected void storeState() {
        if (storePartials) {
            likelihoodCore.storeState();
        }
        super.storeState();
        

    }

    /**
     * Restore the additional stored state
     */
    @Override
	protected void restoreState() {

        if (storePartials) {
            likelihoodCore.restoreState();
        } else {
            updateAllNodes();
        }

        super.restoreState();

    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    @Override
	protected double calculateLogLikelihood() {

        if (patternLogLikelihoods == null) {
            patternLogLikelihoods = new double[patternCount];
        }

        final NodeRef root = treeModel.getRoot();
        traverse(treeModel, root);
        
        double logL = 0.0;
        double ascertainmentCorrection = getAscertainmentCorrection(patternLogLikelihoods);
        for (int i = 0; i < patternCount; i++) {
            logL += (patternLogLikelihoods[i] - ascertainmentCorrection) * patternWeights[i];
        }
        
        if (logL == Double.NEGATIVE_INFINITY) {
            Logger.getLogger("dr.evomodel").info("TreeLikelihood, " + this.getId() + ", turning on partial likelihood scaling to avoid precision loss");

            // We probably had an underflow... turn on scaling
            likelihoodCore.setUseScaling(true);
            
            // and try again...
            updateAllNodes();
            updateAllPatterns();
            traverse(treeModel, root);

            logL = 0.0;
            ascertainmentCorrection = getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < patternCount; i++) {
                logL += (patternLogLikelihoods[i] - ascertainmentCorrection) * patternWeights[i];
            }
        }

        //********************************************************************
        // after traverse all nodes and patterns have been updated --
        //so change flags to reflect this.
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = false;
        }
        //********************************************************************

        return logL;
    }

    public double[] getPatternLogLikelihoods() {
        getLogLikelihood(); // Ensure likelihood is up-to-date
        double ascertainmentCorrection = getAscertainmentCorrection(patternLogLikelihoods);
        double[] out = new double[patternCount];
        for (int i = 0; i < patternCount; i++) {
            if (patternWeights[i] > 0) {
                out[i] = (patternLogLikelihoods[i] - ascertainmentCorrection) * patternWeights[i];
            } else {
                out[i] = Double.NEGATIVE_INFINITY;
            }
        }
        return out;
    }

    /* Calculate ascertainment correction if working off of AscertainedSitePatterns
    @param patternLogProbs log pattern probabilities
    @return the log total probability for a pattern.
    */
    protected double getAscertainmentCorrection(double[] patternLogProbs) {
//        if (patternList instanceof AscertainedSitePatterns) {
//            return ((AscertainedSitePatterns) patternList).getAscertainmentCorrection(patternLogProbs);
//        } else {
            return 0.0;
//        }
    }

    /**
     * Check whether the scaling is still required. If the sum of all the logScalingFactors
     * is zero then we simply turn off the useScaling flag. This will speed up the likelihood
     * calculations when scaling is not required.
     */
    public void checkScaling() {
//	    if (useScaling) {
//	        if (scalingCheckCount % 1000 == 0) {
//	            double totalScalingFactor = 0.0;
//	            for (int i = 0; i < nodeCount; i++) {
//	                for (int j = 0; j < patternCount; j++) {
//	                    totalScalingFactor += scalingFactors[currentPartialsIndices[i]][i][j];
//	                }
//	            }
//	            useScaling = totalScalingFactor < 0.0;
//	            Logger.getLogger("dr.evomodel").info("LikelihoodCore total log scaling factor: " + totalScalingFactor);
//	            if (!useScaling) {
//	                Logger.getLogger("dr.evomodel").info("LikelihoodCore scaling turned off.");
//	            }
//	        }
//	        scalingCheckCount++;
//	    }
    }


    /**
     * Traverse the tree calculating partial likelihoods.
     *
     * @return whether the partials for this node were recalculated.
     */
    protected boolean traverse(Tree tree, NodeRef node) {

        boolean update = false;
        int nodeNum = node.getNumber();
        
        NodeRef parent = tree.getParent(node);
        
        // First update the transition probability matrix(ices) for this branch
        if (parent != null && updateNode[nodeNum]) {

            final double branchRate = branchRateModel.getBranchRate(tree, node);

            // Get the operational time of the branch
            final double branchTime = branchRate * (tree.getNodeHeight(parent) - tree.getNodeHeight(node));

            if (branchTime < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchTime);
            }
            likelihoodCore.setNodeMatrixForUpdate(nodeNum);

            for (int i = 0; i < categoryCount; i++) {

                double branchLength = siteModel.getRateForCategory(i) * branchTime;
                siteModel.getSubstitutionModel().getTransitionProbabilities(branchLength, probabilities);
                likelihoodCore.setNodeMatrix(nodeNum, i, probabilities);
            }

            update = true;
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            final boolean update1 = traverse(tree, child1);

            NodeRef child2 = tree.getChild(node, 1);
            final boolean update2 = traverse(tree, child2);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                final int childNum1 = child1.getNumber();
                final int childNum2 = child2.getNumber();

                likelihoodCore.setNodePartialsForUpdate(nodeNum);
                likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum);

                if (COUNT_TOTAL_OPERATIONS) {
                    totalOperationCount ++;
                }

                if (parent == null) {
                    // No parent this is the root of the tree -
                    // calculate the pattern likelihoods
                    double[] frequencies = frequencyModel.getFrequencies();
                    double[] partials = getRootPartials();

                    likelihoodCore.calculateLogLikelihoods(partials, frequencies, patternLogLikelihoods);
                }
                update = true;
            }
        }
        else{
        	
        	if(nodeNum == updateExternalNodeIndex){
        		update = true;
        	}
        }

        return update;

    }

    public final double[] getRootPartials() {
        if (rootPartials == null) {
            rootPartials = new double[patternCount * stateCount];
        }

        int nodeNum = treeModel.getRoot().getNumber();

        // moved this call to here, because non-integrating siteModels don't need to support it - AD
        double[] proportions = siteModel.getCategoryProportions();
        likelihoodCore.integratePartials(nodeNum, proportions, rootPartials);

        return rootPartials;
    }

    /**
     * the root partial likelihoods (a temporary array that is used
     * to fetch the partials - it should not be examined directly -
     * use getRootPartials() instead).
     */
    private double[] rootPartials = null;
/*
    public class SiteLikelihoodsStatistic extends Statistic.Abstract {

        public SiteLikelihoodsStatistic() {
            super("siteLikelihoods");
        }

        public int getDimension() {
            if (patternList instanceof SitePatterns) {
                return ((SitePatterns)patternList).getSiteCount();
            } else {
                return patternList.getPatternCount();
            }
        }

        public String getDimensionName(int dim) {
            return getTreeModel().getId() + "site-" + dim;
        }

        public double getStatisticValue(int i) {

            if (patternList instanceof SitePatterns) {
                int index = ((SitePatterns)patternList).getPatternIndex(i);
                if( index >= 0 ) {
                    return patternLogLikelihoods[index] / patternWeights[index];
                } else {
                    return 0.0;
                }
            } else {
                return patternList.getPatternCount();
            }
        }

    }

*/    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    /**
     * the frequency model for these sites
     */
    protected final FrequencyModel frequencyModel;

    /**
     * the site model for these sites
     */
    protected final SiteModel siteModel;

    /**
     * the branch rate model
     */
    protected final BranchRateModel branchRateModel;

    /**
     * the tip partials model
     */
//    private final TipStatesModel tipStatesModel;

    private final boolean storePartials;

    protected final boolean integrateAcrossCategories = false;

    /**
     * the categories for each site
     */
    protected int[] siteCategories = null;


    /**
     * the pattern likelihoods
     */
    protected double[] patternLogLikelihoods = null;

    /**
     * the number of rate categories
     */
    protected int categoryCount;

    /**
     * an array used to transfer transition probabilities
     */
    protected double[] probabilities;


    /**
     * an array used to transfer tip partials
     */
//    protected double[] tipPartials;

    /**
     * the LikelihoodCore
     */
    protected LikelihoodCore likelihoodCore;
 
    public String diagnostic(){
    	StringBuilder sb = new StringBuilder();
    	sb.append("Diagnostic!\n");
    	SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int spectrumIndex = record.getSpectrumIndex();
		int siteIndex = record.getSingleIndex();
//    	updateExternalNodeIndex = -1;
//    	spectrumIndex -> taxonName -> indexOnTree
		String taxonId = spectrumModel.getTaxonId(spectrumIndex);
        updateExternalNodeIndex = treeModel.getTaxonIndex(taxonId );
    	
    	
    	for (int i = 0; i < nodeCount; i++) {
    		double[] outPartials = new double[stateCount*patternCount];
    		likelihoodCore.getPartials(i, outPartials);
    		sb.append("Node: "+i +"\t"+ Arrays.toString(outPartials)+"\n");
		}
		
//		sb.append(spectrumIndex +"\t"+ siteIndex +"\tTreeNode:"+ updateExternalNodeIndex +"\n" );
//		sb.append(taxonId +"\t"+ treeModel.getTaxonId(updateExternalNodeIndex) +"\n");
//		spectrumModel.getSpecturmFrequencies(spectrumIndex, siteIndex);
//		sb.append(Arrays.toString(outPartials));
		
		
		sb.append("End diagnostic!\n");
    	return sb.toString();
    }

	public static String compareTwoModels(
			SpectrumTreeLikelihood spectrumTreeLikelihood,
			SpectrumTreeLikelihood newSpectrumTreeLikelihood) {

		StringBuilder sb = new StringBuilder();
		sb.append(Arrays.toString(spectrumTreeLikelihood.getPatternLogLikelihoods()));
		sb.append(Arrays.toString(newSpectrumTreeLikelihood.getPatternLogLikelihoods()));
		
		int nodeCount = spectrumTreeLikelihood.nodeCount;
		int stateCount = 4;
		int patternCount = spectrumTreeLikelihood.patternCount;
		LikelihoodCore likelihoodCore = spectrumTreeLikelihood.getLikelihoodCore();
		LikelihoodCore likelihoodCore2 = newSpectrumTreeLikelihood.getLikelihoodCore();
		for (int i = 0; i < nodeCount; i++) {
    		double[] outPartials = new double[stateCount*patternCount];
    		double[] outPartials2 = new double[stateCount*patternCount];
    		likelihoodCore.getPartials(i, outPartials);
    		likelihoodCore2.getPartials(i, outPartials2);
    		for (int j = 0; j < outPartials2.length; j++) {
				if(outPartials[j]!=outPartials2[j]){
					sb.append("node: "+i+" index: "+j +"\t"+ outPartials[j] +"\t"+ outPartials2[j]+"\n");
				}
			}
//    		sb.append("Node: "+i +"\t"+ Arrays.toString(outPartials)+"\n");
		}
		return sb.toString();
	}
}
