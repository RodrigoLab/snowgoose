package srp.likelihood;


import java.util.ArrayList;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.HaplotypeModel;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequences;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree.MissingTaxonException;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.inference.model.Parameter;

public class LikelihoodCalculation {//extends TraceCorrelationAssert {

	
//	protected static final String TREE_HEIGHT = TreeModel.TREE_MODEL + "." + TreeModelParser.ROOT_HEIGHT;
    protected TreeModel treeModel;
    protected Alignment alignment;
    protected BranchRateModel branchRateModel;
    
    protected AlignmentMapping aMap;
    protected HaplotypeModel haplotypeModel;
    
    @Deprecated
    protected Sequences shortReads;

    private static GammaSiteModel siteModel;
    
    private CoalescentLikelihood coalescent;
    private TreeLikelihood treeLikelihood;
	
    private ShortReadLikelihood shortReadLikelihood;
	private int LogLikelihood;
	
	public LikelihoodCalculation() {
	}

	// public LikelihoodCalculation(FlexibleTree tree, Alignment
	// alignment) {
	// setTreeAndAlignment(tree, alignment);
	//
	// }

	public LikelihoodCalculation(TreeModel treeModel, Parameter popSize) {
		this.treeModel = treeModel;
		setPopSize(popSize);
	}

	public LikelihoodCalculation(TreeModel treeModel, Alignment alignment) {
		setTreeAndAlignment(treeModel, alignment);
	}

	public LikelihoodCalculation(TreeModel treeModel,
			Alignment alignment, Parameter popSize) {
		this(treeModel, alignment);
		setPopSize(popSize);
	}

	public LikelihoodCalculation(TreeModel treeModel,
			Alignment alignment, Parameter popSize, Parameter mutationRate) {
		this(treeModel, alignment);
		setPopSize(popSize);
		setMutationRate(mutationRate);
	}

	@Deprecated
	public LikelihoodCalculation(TreeModel treeModel,
			Alignment alignment, Sequences shortReads) {
		this(treeModel, alignment);
		setShortReads(shortReads);
	}
	
	
	public LikelihoodCalculation(TreeModel treeModel, AlignmentMapping aMap){
		//TODO implement later, auot g
	}
	
	public LikelihoodCalculation(TreeModel treeModel, HaplotypeModel haplotypeModel){
//		Alignment alignment, Alignment shortReads) {
	
		alignment = haplotypeModel.getAlignment();
		setTreeAndAlignment(treeModel, alignment);
	}
	public LikelihoodCalculation(TreeModel treeModel, AlignmentMapping aMap, HaplotypeModel haplotypeModel){
//			Alignment alignment, Alignment shortReads) {
		
		alignment = haplotypeModel.getAlignment();
		setTreeAndAlignment(treeModel, alignment);
		setShortReads(aMap, haplotypeModel);
	}
	
	private void setShortReads(AlignmentMapping aMap, HaplotypeModel haplotypeModel) {
		this.aMap = aMap;
		this.haplotypeModel =haplotypeModel;
		shortReadLikelihood = new ShortReadLikelihood(this.aMap, this.haplotypeModel);
		
	}

	public void setPopSize(double p, double lower, double upper) {
		Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, p, lower, upper);
		setPopSize(popSize);
	}

	public void setPopSize(Parameter popSize) {
		ConstantPopulationModel constantModel = new ConstantPopulationModel(popSize, Units.Type.DAYS);
    	
		try {
			coalescent = new CoalescentLikelihood(treeModel, null, new ArrayList<TaxonList>(), constantModel);
			coalescent.setId("coalescent");
		} catch (MissingTaxonException e) {
			e.printStackTrace();
		}
    	
		
	}

	public void setMutationRate(double u, double lower, double upper){
		Parameter rateParameter =  new Parameter.Default(StrictClockBranchRates.RATE, u, lower, upper);
		setMutationRate(rateParameter);
	}
	
	public void setMutationRate(Parameter u){
    	branchRateModel = new StrictClockBranchRates(u);
    	setupJCTreeLikelihood();
	}
	
	public void setTreeAndAlignment(FlexibleTree tree, Alignment alignment) {
    	TreeModel treeModel = new TreeModel(TreeModel.TREE_MODEL, tree, false, false);
		setTreeAndAlignment(treeModel, alignment);
	}
	
    public void setTreeAndAlignment(TreeModel treeModel, Alignment alignment) {
		this.treeModel = treeModel;
		this.alignment = alignment;
		setupJCTreeLikelihood();
		
		
	}
    
//	private void setupJCTreeLikelihood(BranchRateModel branchRateModel) {
//		
//    	SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
//
//		GammaSiteModel siteModel = setupDefaultJC69SiteModel();
//		treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel,
//				branchRateModel, null, false, false, true, false, false);
//		
//	}
	
    @Deprecated
	public void setShortReads(Sequences shortReads) {
		this.shortReads = shortReads;
		shortReadLikelihood = new ShortReadLikelihood(shortReads, alignment);
		
	}

	private void setupJCTreeLikelihood() {
		
    	SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
    	
		GammaSiteModel siteModel = setupDefaultJC69SiteModel();
		
		treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel,
				branchRateModel, null, false, false, true, false, false);

		
		//		treeLikelihood.modelChangedEvent(model, object, index)
//		System.out.println(treeLikelihood.getModelCount());
//		for (int i = 0; i < treeLikelihood.getModelCount(); i++) {
//			System.out.println("z"+treeLikelihood.getModel(i).getModelName());
//		}
		
	}

	public double getTreeLikelihood(){
    	
    	return treeLikelihood.getLogLikelihood();

    }
    public double getShortReadLikelihood(){
		return shortReadLikelihood.getLogLikelihood();
	}

	public double getCoalescentLikelhood() {
		return coalescent.getLogLikelihood();
		
	}
	
	public double getLoglikelihood(){
		LogLikelihood = 0;
//		LogLikelihood += getTreeLikelihood();
		LogLikelihood += getShortReadLikelihood();
//		LogLikelihood += getCoalescentLikelhood();
		
		return LogLikelihood;
	}
//
//	public double calPopLikelihood() {
//
//    	Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 380.0, 0, 38000.0);
//
//    	////////
//
//    	ConstantPopulationModel startingTree = new ConstantPopulationModel(popSize, Units.Type.YEARS);
//    	//    	        ConstantPopulation constant = (ConstantPopulation) startingTree.getDemographicFunction();
//
//    	//////
//    	//    	        CoalescentSimulator simulator = new CoalescentSimulator();
//    	//    	        Tree tree = simulator.simulateTree(alignment, constant);
//    	//    	        treeModel = new TreeModel(tree);//treeModel
//
//
//    	////////
//    	ConstantPopulationModel constantModel = startingTree;
//
//    	CoalescentLikelihood coalescent = null;
//    	try {
//    		coalescent = new CoalescentLikelihood(treeModel, null, new ArrayList<TaxonList>(), constantModel);
//    		coalescent.setId("coalescent");
//    	} catch (MissingTaxonException e) {
//    		e.printStackTrace();
//    	}
//    	//        System.out.println(coalescent.calculateLogLikelihood());
//    	System.out.println(coalescent.getLogLikelihood());
//    	// clock model
//    	Parameter rateParameter =  new Parameter.Default(StrictClockBranchRates.RATE, 1e-5, 0, 100.0);
//    	StrictClockBranchRates branchRateModel = new StrictClockBranchRates(rateParameter);
//
//    	GammaSiteModel siteModel = setupDefaultJC69SiteModel();
//    	SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
//    	TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, 
//    			branchRateModel, null, false, false, true, false, false);
//         
//		return treeLikelihood.getLogLikelihood();
//    }
    
    
    public static GammaSiteModel setupDefaultJC69SiteModel(){
    	
    	Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25});
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 1, 1);
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        siteModel = new GammaSiteModel(hky);
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteModel.setMutationRateParameter(mu);
    	
        return siteModel;
    }

	public void updateHaplotypes(HaplotypeModel haplotypeModel) {
		this.haplotypeModel = haplotypeModel;
//		alignment = this.haplotypes.getAlignment();
		shortReadLikelihood.updateHaplotypes(this.haplotypeModel);
		
//		patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
//		setTreeAndAlignment(treeModel, alignment);
		
//		SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
//    	
//		GammaSiteModel siteModel = setupDefaultJC69SiteModel();
//		
//		treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel,
//				branchRateModel, null, false, false, true, false, false);

		
	}

	public void storeState(){
		shortReadLikelihood.storeState();
		
	}
	
	public void restorreState(){
		shortReadLikelihood.restoreState();
		
	}

	public void calculateShortReadLikelihoodFull() {
		shortReadLikelihood.calculateLogLikelihoodSelect(0);
		
	}
	
}