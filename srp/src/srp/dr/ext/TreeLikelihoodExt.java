package srp.dr.ext;

import srp.evolution.OperationRecord;
import srp.evolution.haplotypes.old.OldHaplotypeModel;
import srp.haplotypes.HaplotypeModel;
import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evolution.util.TaxonList.MissingTaxonException;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AbstractLikelihoodCore;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.inference.model.Model;

public class TreeLikelihoodExt extends TreeLikelihood {

	private HaplotypeModel haplotypeModel;
	private SitePatternsExt sitePatternExt;
	private int updateExternalNodeIndex;
	private int[] tempstates;
	private AbstractLikelihoodCore likelihoodCoreA;
	
	public TreeLikelihoodExt(HaplotypeModel haplotypeModel, TreeModel treeModel,
			SiteModel siteModel, BranchRateModel branchRateModel,
			TipStatesModel tipStatesModel, boolean useAmbiguities,
			boolean allowMissingTaxa, boolean storePartials,
			boolean forceJavaCore, boolean forceRescaling) {
		
		
		super(new SitePatternsExt (haplotypeModel),
				treeModel, siteModel, branchRateModel,
				tipStatesModel, useAmbiguities, allowMissingTaxa,
				storePartials, forceJavaCore, forceRescaling);
		
		this.sitePatternExt = (SitePatternsExt) getPatternList(); 
		this.haplotypeModel = haplotypeModel;
		addModel(this.haplotypeModel);
		tempstates = new int[patternCount];
		likelihoodCoreA = (AbstractLikelihoodCore)likelihoodCore;
//		LogManager.getLogManager().reset();
//		final Logger logger = Logger.getLogger("dr.evomodel");
//        logger.setLevel(Level.OFF);
        
	}

	/**
	 * 
	 */
	private static final long serialVersionUID = 6348244745369755426L;

    public static long time1 = 0;
    public static long time2 = 0;
    public static long time3 = 0;
    /**
     * Handles model changed events from the submodels.
     */
    @Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
//    	System.out.println("handleModelChangedEvent in TreeLikelihoodExt\t"+model.getModelName());
    	if (model == haplotypeModel){

    		sitePatternExt.updateAlignment(haplotypeModel);
    		updatePatternListExt();
    		likelihoodKnown = false;
    	}
    	else if (model == oldHaplotypeModel){ //REMOVE: Remove OldHaplotype
    		System.out.println("BAD! using oldHaplotypeModel");
//    		sitePatternExt.updateAlignment(oldHaplotypeModel);
//    		updatePatternListExt(sitePatternExt);
    		likelihoodKnown = false;
    	}
    	else{
    		super.handleModelChangedEvent(model, object, index);
    	}
    }
    @Override
	protected void restoreState() {
//    	System.out.println("restore state");
//        if (storePartials) {
//            likelihoodCore.restoreState();
//        } else {
//            updateAllNodes();
//        }
//    	sitePatternExt.updateAlignment(haplotypeModel);
//		updatePatternListExt(sitePatternExt);
        super.restoreState();

    }
    
    @Override
    protected void storeState(){
//    	System.out.println("storeState");
    	super.storeState();
    }


	public void updatePatternListExt() {

//        this.patternList = sitePatternExt;
//patternList
        OperationRecord record = haplotypeModel.getOperationRecord();
		int haplotypeIndex = record.getSpectrumIndex();
        String taxonId = haplotypeModel.getTaxonId(haplotypeIndex);
		updateExternalNodeIndex = treeModel.getTaxonIndex(taxonId);

		updateNode[updateExternalNodeIndex] = true;
		int site;
//        for (int i = 0; i < nodeCount; i++) {
//            updateNode[i] = true;
//        }
		likelihoodCoreA.getNodeStates(updateExternalNodeIndex, tempstates);

		switch (record.getOperation()) {
		case SINGLE:
			site = record.getSingleIndex();
			tempstates[site] = patternList.getPatternState(haplotypeIndex, site);
			likelihoodCore.setNodeStates(updateExternalNodeIndex, tempstates);
			break;
		case MULTI:
			int[] sites = record.getAllSiteIndexs();
			for (int s : sites) {
				tempstates[s] = patternList.getPatternState(haplotypeIndex, s);
			}
			likelihoodCore.setNodeStates(updateExternalNodeIndex, tempstates);
			break;
		case COLUMN:
			site = record.getSingleIndex();
			for (int h = 0; h < haplotypeModel.getHaplotypeCount(); h++) {
				haplotypeIndex = h;
		        taxonId = haplotypeModel.getTaxonId(haplotypeIndex);
				updateExternalNodeIndex = treeModel.getTaxonIndex(taxonId);

				updateNode[updateExternalNodeIndex] = true;
				likelihoodCoreA.getNodeStates(updateExternalNodeIndex, tempstates);
				tempstates[site] = patternList.getPatternState(haplotypeIndex, site);
				likelihoodCore.setNodeStates(updateExternalNodeIndex, tempstates);

			}
			
			
			break;
		case RECOMBINATION:
			int[] twoPositions = record.getRecombinationPositionIndex();
			int[] twoHapIndexs = record.getRecombinationSpectrumIndex();
			for (int h : twoHapIndexs) {
				haplotypeIndex = h;
		        taxonId = haplotypeModel.getTaxonId(haplotypeIndex);
				updateExternalNodeIndex = treeModel.getTaxonIndex(taxonId);

				updateNode[updateExternalNodeIndex] = true;
				likelihoodCoreA.getNodeStates(updateExternalNodeIndex, tempstates);
			
				for (int s = twoPositions[0]; s < twoPositions[1]; s++) {
					tempstates[s] = patternList.getPatternState(haplotypeIndex, s);
				}
				likelihoodCore.setNodeStates(updateExternalNodeIndex, tempstates);

			}
			break;
		default:
//			break;
			throw new IllegalArgumentException("Invalid operation type:"+record.getOperation());
		}
//		int[] sites = record.getAllSiteIndexs();
//		((AbstractLikelihoodCore)likelihoodCore).getNodeStates(updateExternalNodeIndex, tempstates);
//		for (int s : sites) {
//			tempstates[s] = patternList.getPatternState(haplotypeIndex, s);
//		}
		

	}

//	public void updatePatternList(SitePatternsExt patterns, OldHaplotypeModel haplotypeModel) {
////		TODO: more test required
////		Alignment alignment = haplotypeModel.getAlignment();
//        
//		patterns.updateAlignment(haplotypeModel);
//		updatePatternListExt(patterns);
//	}
//
//	public void updatePatternList(OldHaplotypeModel haplotypeModel) {
////		TODO: more test required
////		sitePatternExt.updateAlignment(haplotypeModel);
//		updatePatternListExt(haplotypeModel);
//	}

	public SiteModel getSiteModel(){
		return siteModel;
	}
	

	@Override
	protected boolean traverse(Tree tree, NodeRef node) {

        boolean update = false;

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);
//System.out.print(nodeNum);
        // First update the transition probability matrix(ices) for this branch
        if (parent != null && updateNode[nodeNum]) {
//System.out.print("\tupdate trans prob\t");
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
//        	System.out.print("\tupdate internal\t");
            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            final boolean update1 = traverse(tree, child1);

            NodeRef child2 = tree.getChild(node, 1);
            final boolean update2 = traverse(tree, child2);

            // If either child node was updated then update this node too
            if (update1 || update2) {
//            	System.out.print("\t update partial\t");
                final int childNum1 = child1.getNumber();
                final int childNum2 = child2.getNumber();

                likelihoodCore.setNodePartialsForUpdate(nodeNum);

                if (integrateAcrossCategories) {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum);
                } else {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (COUNT_TOTAL_OPERATIONS) {
                    totalOperationCount ++;
                }

                if (parent == null) {
                    // No parent this is the root of the tree -
                    // calculate the pattern likelihoods
                    double[] frequencies = frequencyModel.getFrequencies();

                    double[] partials = getRootPartialsExt(); //TODO deal with private final

                    likelihoodCore.calculateLogLikelihoods(partials, frequencies, patternLogLikelihoods);
                }

                update = true;
            }
        }
 else{
        	
        	if(nodeNum == updateExternalNodeIndex){
//        		System.out.println("update:"+updateExternalNodeIndex);
        		update = true;
        	}
        }
//System.out.println();
        return update;

    }

    public final double[] getRootPartialsExt() {
        if (rootPartials == null) {
            rootPartials = new double[patternCount * stateCount];
        }

        int nodeNum = treeModel.getRoot().getNumber();
        if (integrateAcrossCategories) {

            // moved this call to here, because non-integrating siteModels don't need to support it - AD
            double[] proportions = siteModel.getCategoryProportions();
            likelihoodCore.integratePartials(nodeNum, proportions, rootPartials);
        } else {
            likelihoodCore.getPartials(nodeNum, rootPartials);
        }

        return rootPartials;
    }

    /**
     * the root partial likelihoods (a temporary array that is used
     * to fetch the partials - it should not be examined directly -
     * use getRootPartials() instead).
     */
    private double[] rootPartials = null;
	
    @Deprecated
	private void resetRootPartials() {
		rootPartials = null;		
	}

	@Deprecated private OldHaplotypeModel oldHaplotypeModel;
	
	@Deprecated
		public void updatePatternListExt(PatternList patternList) {
			
	
	        this.patternList = patternList;
	//        this.dataType = patternList.getDataType();
	        patternCount = patternList.getPatternCount();
	//        stateCount = dataType.getStateCount();
	
	        patternWeights = patternList.getPatternWeights();
	
	//        this.treeModel = treeModel;
	//        addModel(treeModel);
	
	//        nodeCount = treeModel.getNodeCount();
	
	//        updateNode = new boolean[nodeCount];
	        
	        for (int i = 0; i < nodeCount; i++) {
	            updateNode[i] = true;
	        }
	
	        likelihoodKnown = false;
			
	        patternLogLikelihoods = new double[patternCount];
	        resetRootPartials(); 
	//        getRootPartials();
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
	//		boolean forceRescaling = false;
	
			        try {
	//		            this.siteModel = siteModel;
	//		            addModel(siteModel);
	
	//		            this.frequencyModel = siteModel.getFrequencyModel();
	//		            addModel(frequencyModel);
	
	//		            this.tipStatesModel = tipStatesModel;
	
	//		            integrateAcrossCategories = siteModel.integrateAcrossCategories();
	
			            categoryCount = siteModel.getCategoryCount();
	
	//		            final Logger logger = Logger.getLogger("dr.evomodel");
	//		            String coreName = "Java general";
	//		            boolean forceJavaCore = false;
	//		            if (integrateAcrossCategories) {
	//
	//		                final DataType dataType = patternList.getDataType();
	//
	//		                if (dataType instanceof dr.evolution.datatype.Nucleotides) {
	//
	//		                    if (!forceJavaCore && NativeNucleotideLikelihoodCore.isAvailable()) {
	//		                        coreName = "native nucleotide";
	//		                        likelihoodCore = new NativeNucleotideLikelihoodCore();
	//		                    } else {
	//		                        coreName = "Java nucleotide";
	//		                        likelihoodCore = new NucleotideLikelihoodCore();
	//		                    }
	//
	//		                } 
	//		            } else {
	//		                likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
	//		            }
	
			            probabilities = new double[stateCount * stateCount];
	
			            likelihoodCore.initialize(nodeCount, patternCount, categoryCount, integrateAcrossCategories);
	
			            int extNodeCount = treeModel.getExternalNodeCount();
			            int intNodeCount = treeModel.getInternalNodeCount();
			           
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
	
	
	
			        } catch (TaxonList.MissingTaxonException mte) {
			            throw new RuntimeException(mte.toString());
			        }
	
	//		        addStatistic(new SiteLikelihoodsStatistic());
	//		        System.out.println(getStatisticCount());
	//		        System.out.println(getStatistic(0).getDimension()+"\t"+ getStatistic(0).getStatisticValue(10));
		}
	@Deprecated
	public TreeLikelihoodExt(OldHaplotypeModel haplotypeModel, TreeModel treeModel,
			SiteModel siteModel, BranchRateModel branchRateModel,
			TipStatesModel tipStatesModel, boolean useAmbiguities,
			boolean allowMissingTaxa, boolean storePartials,
			boolean forceJavaCore, boolean forceRescaling) {
		
		
		super(new SitePatternsExt (haplotypeModel, null, 0, -1, 1, true),
				treeModel, siteModel, branchRateModel,
				tipStatesModel, useAmbiguities, allowMissingTaxa,
				storePartials, forceJavaCore, forceRescaling);
		
		this.sitePatternExt = (SitePatternsExt) getPatternList(); 
		this.oldHaplotypeModel = haplotypeModel;
		addModel(this.oldHaplotypeModel);
	}
}
