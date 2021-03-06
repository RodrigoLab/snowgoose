package dr.evomodel.antigenic;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.MCMCOperator;
import dr.math.GammaFunction;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Gabriela Cybis
 * @author Marc Suchard
 */
public class NPAntigenicLikelihood extends AbstractModelLikelihood {
	
	  public static final String NP_ANTIGENIC_LIKELIHOOD = "NPAntigenicLikelihood";

	
	
	public NPAntigenicLikelihood (TreeModel treeModel,CompoundParameter traitParameter , Parameter assignments, Parameter links, Parameter clusterPrec, Parameter priorMean, Parameter priorPrec, double transformFactor,Parameter means1, Parameter means2 ){
		  super(NP_ANTIGENIC_LIKELIHOOD);


          this.assignments = assignments;
          this.links = links;
		  this.clusterPrec = clusterPrec;
		  this.priorPrec = priorPrec;
		  this.priorMean = priorMean;
	      this.treeModel= treeModel;
		  this.traitParameter= traitParameter;
		  this.transformFactor=transformFactor;
		  this.means1=means1;
          this.means2=means2;

		  this.alpha= 1.0;

         addVariable(clusterPrec);

        addVariable(traitParameter);
        addVariable(assignments);
        addVariable(links);
        addModel(treeModel);

        numdata = traitParameter.getParameterCount();
        this.allTips=Tree.Utils.getExternalNodes(treeModel,treeModel.getRoot());


        setData();
       setDepMatrix();



	 for (int i=0; i<numdata; i++){
		 assignments.setParameterValue(i, i);
	     links.setParameterValue(i,i);
     }
	 
	
	       
	       this.logLikelihoodsVector = new double[links.getDimension()+1];

        double[] m = new double[2];
        m[0]= priorMean.getParameterValue(0);
        m[1]= priorMean.getParameterValue(1);


         double v0 = 2;
         double v1 = 3;

        double k0= priorPrec.getParameterValue(0)/clusterPrec.getParameterValue(0);
        double k1= k0+1;


        double[][] T0Inv= new double[2][2];
        T0Inv[0][0]= v0/clusterPrec.getParameterValue(0);
        T0Inv[1][1]= v0/clusterPrec.getParameterValue(0);
        T0Inv[1][0]= 0.0;
        T0Inv[0][1]= 0.0;


          double logDetT0= -Math.log(T0Inv[0][0]*T0Inv[1][1]);


        for(int i=0;i<logLikelihoodsVector.length-1;i++){

         double[][] T1Inv = new double[2][2];
            T1Inv[0][0]=T0Inv[0][0]+(k0/k1)* data[i][0]*data[i][0];
            T1Inv[0][1]=T0Inv[0][1]+(k0/k1)* data[i][0]*data[i][1];
            T1Inv[1][0]=T0Inv[1][0]+(k0/k1)* data[i][1]*data[i][0];
            T1Inv[1][1]=T0Inv[1][1]+(k0/k1)* data[i][1]*data[i][1];


            double logDetT1=-Math.log(T1Inv[0][0]*T1Inv[1][1]-T1Inv[0][1]*T1Inv[1][0]);


            logLikelihoodsVector[i]= -(1*2/2)*Math.log(Math.PI);
            logLikelihoodsVector[i]+= Math.log(k0) - Math.log(k1);
            logLikelihoodsVector[i]+= (v1/2)*logDetT1 - (v0/2)*logDetT0;
            logLikelihoodsVector[i]+= GammaFunction.lnGamma(v1/2)+ GammaFunction.lnGamma((v1/2)-0.5);
            logLikelihoodsVector[i]+=-GammaFunction.lnGamma(v0/2)- GammaFunction.lnGamma((v0/2)-0.5);


        }


       printInformtion(logLikelihoodsVector[0]);
        printInformtion(logLikelihoodsVector[2]);




    }



    private void setData(){
        int dim = traitParameter.getParameter(0).getSize();


        double Data[][] =new double[numdata][dim];

        for (int i=0; i<numdata; i++){
            for (int j=0; j<dim; j++){
                Data[i][j]= traitParameter.getParameter(i).getParameterValue(j);
            }
        }

        this.data=Data;
    }
	

    private void setDepMatrix(){

        depMatrix=new double[numdata][numdata];
        List<NodeRef> childList = new ArrayList<NodeRef>();

        recursion(treeModel.getRoot(),childList);
        logCorrectMatrix(transformFactor);
        // printInformtion(depMatrix);
        logDepMatrix =  new double[numdata][numdata];
        for(int i=0;i<numdata;i++){
            for(int j=0;j<i;j++){
                logDepMatrix[i][j]=Math.log(depMatrix[i][j]);
                logDepMatrix[j][i]=logDepMatrix[j][i];

            }
        }


    }




	
	
	  public Model getModel() {
	        return this;
	    }

	  
	  
	  public double[] getLogLikelihoodsVector(){
		  return logLikelihoodsVector;
	  }

    public Parameter getLinks(){
        return links;
    }

    public Parameter getAssignments(){
        return assignments;
    }


    public double[][] getData(){
		  return data;
	  }
	  
	  public double[][] getDepMatrix(){
		  return depMatrix;
	  }
	  

	  public double[][] getLogDepMatrix(){
		  return logDepMatrix;
	  }
	  
	  public Parameter getPriorMean(){
		  return priorMean;
	  }
	  public Parameter getPriorPrec(){
		  return priorPrec;
	  }
	  public Parameter getClusterPrec(){
		  return clusterPrec;
	  }
	  
	  public void setLogLikelihoodsVector(int pos, double value){
		  logLikelihoodsVector[pos]=value;
	  }

    public void setAssingments(int pos, double value){
        assignments.setParameterValue(pos,value);
    }

    public void setLinks(int pos, double value){
        links.setParameterValue(pos,value);
    }

    public void setMeans(int pos, double[] value){
        means1.setParameterValue(pos,value[0]);
        means2.setParameterValue(pos,value[1]);
    }

	  
	  public double getLogLikelihood() {
          setDepMatrix();
          setData();
         //printInformtion(data[0][0]);
          //printInformtion(depMatrix[0][1]);


		  double logL = 0.0;
		  for (int j=0 ; j<logLikelihoodsVector.length;j++){
			  if(logLikelihoodsVector[j]!=0){
				  logL +=logLikelihoodsVector[j];
			  }
		  }
		  
		  for (int j=0 ; j<links.getDimension();j++){
		if(links.getParameterValue(j)==j){
			logL += Math.log(alpha);
		}
		else{logL += Math.log(depMatrix[j][(int) links.getParameterValue(j)]);
			
		}
		
		double sumDist=0.0;
		for (int i=0;i<numdata;i++){
			if(i!=j){sumDist += depMatrix[i][j];
			}
			}
		
		  logL-= Math.log(alpha+sumDist);
		  }
		
	
		  
		 
		 return logL;
	  }
	 
	  
	  
	  /* Marc's suggestion on recursion for getting matrix from tree*/
	  
	  void recursion( NodeRef node, List childList){
		 
		  
		 List<NodeRef> leftChildTipList = new ArrayList<NodeRef>(); 
		 List<NodeRef> rightChildTipList = new ArrayList<NodeRef>(); 
		  
		 if(!treeModel.isExternal(node)){
			 recursion(treeModel.getChild(node, 0),leftChildTipList);
			 recursion(treeModel.getChild(node, 1),rightChildTipList);
			 
			
			 double lBranch = treeModel.getBranchLength(treeModel.getChild(node, 0));
			 double rBranch = treeModel.getBranchLength(treeModel.getChild(node, 1));
			 

			 Set<NodeRef> notLeftChildList = new HashSet<NodeRef>();
			 notLeftChildList.addAll(allTips);
			 	for (NodeRef i :leftChildTipList){
			 			notLeftChildList.remove(i);			 
}
			 	Set<NodeRef> notRightChildList = new HashSet<NodeRef>();
				 notRightChildList.addAll(allTips);
			 	for (NodeRef i :rightChildTipList){
			 		notRightChildList.remove(i);			 
}
			 
			 for (NodeRef lChild : leftChildTipList){
				  for (NodeRef Child : notLeftChildList){
				depMatrix[Child.getNumber()][lChild.getNumber()] += lBranch;	 
				depMatrix[lChild.getNumber()][Child.getNumber()] += lBranch;	 
				 }
			 }


			 for (NodeRef rChild : rightChildTipList){
					 for (NodeRef Child : notRightChildList){
					depMatrix[Child.getNumber()][rChild.getNumber()] += rBranch;	 
					depMatrix[rChild.getNumber()][Child.getNumber()] += rBranch;	 
					 }
				 }
			 
			 
			 childList.addAll(leftChildTipList);
			 childList.addAll(rightChildTipList);
		 }
		 else{
			 childList.add(node);
		 }
		 
	  }
	  
	  
	  
	  void logCorrectMatrix(double p){
		  for (int i=0; i<numdata; i++){
		    	for (int j=0; j<i; j++){
		    		depMatrix[i][j]=1/Math.pow(depMatrix[i][j],p);	    
		    		depMatrix[j][i]=depMatrix[i][j];	    
				    }}
		    
	  }
	  
	  
	  
	  
	  // Slow method for computing matrix from tree
	  
	  
	  
	  public double[][] getMatrixFromTree(double p){
		  double[][] Mat = new double[numdata][numdata];
		  
	  for (int i = 0 ; i<numdata; i++){
		  for (int j =0 ; j<i; j++){
			  Mat[i][j] = -p*Math.log(getTreeDist(i,j));			
			  Mat[j][i] = Mat[i][j];
		  }
	  }
	  return Mat;
	  }
	  
	  
	  
	  
	  public double getTreeDist(int i, int j){
		 double dist=0; 
  
		 NodeRef MRCA = findMRCA(i,j);

		 NodeRef Parent = treeModel.getExternalNode(i);
		 while (Parent!=MRCA){
			 dist+=treeModel.getBranchLength(Parent);
			 Parent = treeModel.getParent(Parent);
		 }

		 Parent = treeModel.getExternalNode(j);
		 while (Parent!=MRCA){
			 dist+=treeModel.getBranchLength(Parent);
			 Parent = treeModel.getParent(Parent);
		 }
		  
			 return dist;
	  }
	  
	  
	  private NodeRef findMRCA(int iTip, int jTip) {
	        Set<String> leafNames = new HashSet<String>();
	        leafNames.add(treeModel.getTaxonId(iTip));
	        leafNames.add(treeModel.getTaxonId(jTip));
	        return Tree.Utils.getCommonAncestorNode(treeModel, leafNames);
	    }


	  public void printInformtion(double[][] Mat) {
         StringBuffer sb = new StringBuffer("matrix \n");
         for(int i=0;i <numdata; i++){
        	 sb.append(" \n");
        	 for(int j=0; j<numdata; j++){
        		 sb.append(Mat[i][j]+" \t");
        	 }
         }
         
         Logger.getLogger("dr.evomodel").info(sb.toString()); };
        
        

   	  public void printOrder() {
            StringBuffer sb = new StringBuffer("taxa \n");
            for(int i=0;i <numdata; i++){
           	 sb.append(" \n");
           	 
           		 sb.append(treeModel.getTaxonId(i));
           	             }
            
            Logger.getLogger("dr.evomodel").info(sb.toString()); };
           
         
         
         
         
         public void printInformtion(double x) {
             StringBuffer sb = new StringBuffer("Info \n");
             		 sb.append(x);
             
             Logger.getLogger("dr.evomodel").info(sb.toString()); };
	  
	
	  
	  
	  
	 
	  
	  
	  
	  public void makeDirty() {
	    }

	    public void acceptState() {
	        // DO NOTHING
	    }

	    public void restoreState() {
	        // DO NOTHING
	    }

	    public void storeState() {
	        // DO NOTHING
	    }

	    protected void handleModelChangedEvent(Model model, Object object, int index) {
	        // DO NOTHING
	    }

	    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
	        // DO NOTHING
	    }

	Set<NodeRef> allTips;  
	CompoundParameter traitParameter;  
	double alpha;
	Parameter clusterPrec ;
	Parameter priorPrec ;
	Parameter priorMean ;
	Parameter assignments;
    Parameter links;
    Parameter means2;
    Parameter means1;
    TreeModel treeModel;
	String traitName;
	double[][] data;
	double[][] depMatrix;
	double[][] logDepMatrix;
	double[] logLikelihoodsVector;
	int numdata;
	double transformFactor;
	
	
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public final static String CLUSTER_PREC = "clusterPrec";
        public final static String PRIOR_PREC = "priorPrec";
        public final static String PRIOR_MEAN = "priorMean";
        public final static String ASSIGNMENTS = "assignments";
        public final static String LINKS = "links";
        public final static String MEANS_1 = "clusterMeans1";
        public final static String MEANS_2 = "clusterMeans2";
        public final static String TRANSFORM_FACTOR = "transformFactor";
        boolean integrate = false;
    	
        
        public String getParserName() {
            return NP_ANTIGENIC_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
     
        	TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        	//String traitName = (String) xo.getAttribute(TRAIT_NAME);
        	
	        XMLObject cxo = xo.getChild(CLUSTER_PREC);
	        Parameter clusterPrec = (Parameter) cxo.getChild(Parameter.class);

	        cxo = xo.getChild(PRIOR_PREC);
	        Parameter priorPrec = (Parameter) cxo.getChild(Parameter.class);

	        cxo = xo.getChild(PRIOR_MEAN);
	        Parameter priorMean = (Parameter) cxo.getChild(Parameter.class);

	        cxo = xo.getChild(ASSIGNMENTS);
	        Parameter assignments = (Parameter) cxo.getChild(Parameter.class);


            cxo = xo.getChild(LINKS);
            Parameter links = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(MEANS_2);
            Parameter means2 = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(MEANS_1);
            Parameter means1 = (Parameter) cxo.getChild(Parameter.class);




            double transformFactor=1.0;
	        if(xo.hasAttribute(TRANSFORM_FACTOR)){
	        	transformFactor = xo.getDoubleAttribute(TRANSFORM_FACTOR);
	        }
	     
	        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();
            String traitName = TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;
         
            
            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                    utilities.parseTraitsFromTaxonAttributes(xo, traitName, treeModel, integrate);
           // traitName = returnValue.traitName;
            CompoundParameter traitParameter = returnValue.traitParameter;
            
	     
	        
	        
	        
	       
	        
	        return new NPAntigenicLikelihood(treeModel,traitParameter,  assignments, links, clusterPrec, priorMean,priorPrec,transformFactor, means1,means2);
	    }

	    //************************************************************************
	    // AbstractXMLObjectParser implementation
	    //************************************************************************

	    public String getParserDescription() {
	        return "conditional likelihood ddCRP";
	    }

	    public Class getReturnType() {
	        return NPAntigenicLikelihood.class;
	    }

	    public XMLSyntaxRule[] getSyntaxRules() {
	        return rules;
	    }

	    private final XMLSyntaxRule[] rules = {
	    		 new StringAttributeRule(TreeTraitParserUtilities.TRAIT_NAME, "The name of the trait for which a likelihood should be calculated"),
	    	
	    		AttributeRule.newDoubleRule(TRANSFORM_FACTOR,true,"p in transformation of distances -p*log(dist)"),
	    		 
	    		 new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
	                        new ElementRule(Parameter.class)
	                }),   
	    		 new ElementRule(PRIOR_PREC,
		    	                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
		           new ElementRule(CLUSTER_PREC,
			                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
   			       new ElementRule(PRIOR_MEAN,
		    	                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(ASSIGNMENTS,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(LINKS,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(MEANS_1,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(MEANS_2,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
				    	                    new ElementRule(TreeModel.class),
				    	                    
	    };
    };


    
  

    String Atribute = null;
	
}
