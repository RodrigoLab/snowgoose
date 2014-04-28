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

package srp.tdmcmc.likelihood;

import java.util.logging.Logger;

import srp.tdmcmc.SuperBeastTreeLikelihood;
import dr.evolution.alignment.AscertainedSitePatterns;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AminoAcidLikelihoodCore;
import dr.evomodel.treelikelihood.GeneralLikelihoodCore;
import dr.evomodel.treelikelihood.LikelihoodCore;
import dr.evomodel.treelikelihood.NativeAminoAcidLikelihoodCore;
import dr.evomodel.treelikelihood.NativeGeneralLikelihoodCore;
import dr.evomodel.treelikelihood.NativeNucleotideLikelihoodCore;
import dr.evomodel.treelikelihood.NucleotideLikelihoodCore;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inference.model.Model;
import dr.inference.model.Statistic;

/**
 * TreeLikelihoodModel - implements a Likelihood Function for sequences on a tree.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeLikelihood.java,v 1.31 2006/08/30 16:02:42 rambaut Exp $
 */

public class TransdimensionalTreeLikelihood extends SuperBeastTreeLikelihood {

	public TransdimensionalTreeLikelihood(PatternList patternList,
			TreeModel treeModel, SiteModel siteModel,
			BranchRateModel branchRateModel, TipStatesModel tipStatesModel,
			boolean useAmbiguities, boolean allowMissingTaxa,
			boolean storePartials, boolean forceJavaCore, boolean forceRescaling) {
		super(patternList, treeModel, siteModel, branchRateModel, tipStatesModel,
				useAmbiguities, allowMissingTaxa, storePartials, forceJavaCore,
				forceRescaling);
		// TODO Auto-generated constructor stub
	}

	/**
	 * 
	 */
	private static final long serialVersionUID = 73649663825074580L;
    
}
