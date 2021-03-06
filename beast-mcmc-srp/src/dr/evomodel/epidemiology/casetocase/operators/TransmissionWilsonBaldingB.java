package dr.evomodel.epidemiology.casetocase.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.epidemiology.casetocase.BranchMapModel;
import dr.evomodel.epidemiology.casetocase.CaseToCaseTreeLikelihood;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * Implements the Wilson-Balding branch swapping move if it moves an entire subtree of the transmission tree.
 *
 * @author Matthew Hall
 */

public class TransmissionWilsonBaldingB extends AbstractTreeOperator {

    private final CaseToCaseTreeLikelihood c2cLikelihood;
    public static final String TRANSMISSION_WILSON_BALDING_B = "transmissionWilsonBaldingB";
    private double logq;
    private static final boolean DEBUG = false;
    private final int tipCount;

    public TransmissionWilsonBaldingB(CaseToCaseTreeLikelihood c2cLikelihood, double weight) {
        this.c2cLikelihood = c2cLikelihood;
        setWeight(weight);
        tipCount = c2cLikelihood.getTreeModel().getExternalNodeCount();
    }

    public double doOperation() throws OperatorFailedException {

        if(DEBUG){
            c2cLikelihood.debugOutputTree("BeforeTWWB.nex", false);
        }

        proposeTree();

        if (c2cLikelihood.getTreeModel().getExternalNodeCount() != tipCount) {
            int newCount = c2cLikelihood.getTreeModel().getExternalNodeCount();
            throw new RuntimeException("Lost some tips in modified SPR! (" +
                    tipCount + "-> " + newCount + ")");
        }

        if(DEBUG){
            c2cLikelihood.debugOutputTree("AfterTWWB.nex", false);
        }

        return logq;
    }

    public void proposeTree() throws OperatorFailedException {
        TreeModel tree = c2cLikelihood.getTreeModel();
        BranchMapModel branchMap = c2cLikelihood.getBranchMap();
        NodeRef i;
        double oldMinAge, newMinAge, newRange, oldRange, newAge, q;
        // choose a random node avoiding root, and nodes that are ineligible for this move because they have nowhere to
        // go
        final int nodeCount = tree.getNodeCount();
        do {
            i = tree.getNode(MathUtils.nextInt(nodeCount));
        } while (tree.getRoot() == i || !eligibleForMove(i, tree, branchMap));
        final NodeRef iP = tree.getParent(i);

        //this one can go anywhere

        NodeRef j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
        NodeRef k = tree.getParent(j);

        while ((k != null && tree.getNodeHeight(k) <= tree.getNodeHeight(i)) || (i == j)) {
            j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            k = tree.getParent(j);
        }

        if (iP == tree.getRoot() || j == tree.getRoot()) {
            throw new OperatorFailedException("Root changes not allowed!");
        }

        if (k == iP || j == iP || k == i) throw new OperatorFailedException("move failed");

        final NodeRef CiP = getOtherChild(tree, iP, i);
        NodeRef PiP = tree.getParent(iP);

        newMinAge = Math.max(tree.getNodeHeight(i), tree.getNodeHeight(j));
        newRange = tree.getNodeHeight(k) - newMinAge;
        newAge = newMinAge + (MathUtils.nextDouble() * newRange);
        oldMinAge = Math.max(tree.getNodeHeight(i), tree.getNodeHeight(CiP));
        oldRange = tree.getNodeHeight(PiP) - oldMinAge;
        q = newRange / Math.abs(oldRange);

        // need to account for the random repainting of iP

        if(branchMap.get(PiP.getNumber())!=branchMap.get(CiP.getNumber())){
            q *= 0.5;
        }

        if(branchMap.get(k.getNumber())!=branchMap.get(j.getNumber())){
            q *= 2;
        }

        tree.beginTreeEdit();


        if (j == tree.getRoot()) {

            // 1. remove edges <iP, CiP>
            tree.removeChild(iP, CiP);
            tree.removeChild(PiP, iP);

            // 2. add edges <k, iP>, <iP, j>, <PiP, CiP>
            tree.addChild(iP, j);
            tree.addChild(PiP, CiP);

            // iP is the new root
            tree.setRoot(iP);

        } else if (iP == tree.getRoot()) {

            // 1. remove edges <k, j>, <iP, CiP>, <PiP, iP>
            tree.removeChild(k, j);
            tree.removeChild(iP, CiP);

            // 2. add edges <k, iP>, <iP, j>, <PiP, CiP>
            tree.addChild(iP, j);
            tree.addChild(k, iP);

            //CiP is the new root
            tree.setRoot(CiP);

        } else {
            // 1. remove edges <k, j>, <iP, CiP>, <PiP, iP>
            tree.removeChild(k, j);
            tree.removeChild(iP, CiP);
            tree.removeChild(PiP, iP);

            // 2. add edges <k, iP>, <iP, j>, <PiP, CiP>
            tree.addChild(iP, j);
            tree.addChild(k, iP);
            tree.addChild(PiP, CiP);
        }

        tree.setNodeHeight(iP, newAge);

        tree.endTreeEdit();

        //
        logq = Math.log(q);

        // repaint the parent to match either its new parent or its new child (50% chance of each).

        if(MathUtils.nextInt(2)==0){
            branchMap.set(iP.getNumber(), branchMap.get(k.getNumber()), true);
        } else {
            branchMap.set(iP.getNumber(), branchMap.get(j.getNumber()), true);
        }

        if(DEBUG){
            c2cLikelihood.checkPartitions();
        }

    }

    public String getPerformanceSuggestion() {
        return "Not implemented";
    }

    private boolean eligibleForMove(NodeRef node, TreeModel tree, BranchMapModel branchMap){
        // to be eligible for this move, the node's parent and grandparent, or parent and other child,
        // must be in the same partition (so removing the parent has no effect on the remaining links of the TT),
        // and the node and its parent must be in different partitions (such that the move does not disconnect anything)

        return ((tree.getParent(tree.getParent(node))!=null
                && branchMap.get(tree.getParent(node).getNumber())
                ==branchMap.get(tree.getParent(tree.getParent(node)).getNumber()))
                || branchMap.get(tree.getParent(node).getNumber())==branchMap.get(getOtherChild(tree,
                tree.getParent(node), node).getNumber()))
                && branchMap.get(tree.getParent(node).getNumber())!=branchMap.get(node.getNumber());
    }

    @Override
    public String getOperatorName() {
        return TRANSMISSION_WILSON_BALDING_B + " (" + c2cLikelihood.getTreeModel().getId() +")";
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TRANSMISSION_WILSON_BALDING_B;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final CaseToCaseTreeLikelihood c2cL
                    = (CaseToCaseTreeLikelihood) xo.getChild(CaseToCaseTreeLikelihood.class);
            final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

            return new TransmissionWilsonBaldingB(c2cL, weight);
        }

        // ************************************************************************
        // AbstractXMLObjectParser implementation
        // ************************************************************************

        public String getParserDescription(){
            return "This element represents a Wilson-Balding move operator, such that the transplantation of the " +
                    "phylogenetic subtree is also transplantation of a transmission subtree.";
        }

        public Class getReturnType(){
            return TransmissionWilsonBaldingB.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules;{
            rules = new XMLSyntaxRule[]{
                    AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                    new ElementRule(CaseToCaseTreeLikelihood.class)
            };
        }
    };



}
