package srp.core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import srp.dr.ext.TreeLikelihoodExt;
import srp.haplotypes.HaplotypeModel;
import srp.operator.haplotypes.BaseSingleOperator;
import srp.operator.haplotypes.BasesMultiOperator;
import srp.operator.haplotypes.ColumnOperator;
import srp.operator.haplotypes.HaplotypeSwapSectionOperator;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.DeltaExchangeOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.Scalable;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.UniformOperator;
import dr.inference.operators.UpDownOperator;

public class MCMCSetupHelperHaplotype extends MCMCSetupHelper {

	public static HashMap<String, Object> setupTreeLikelihoodHaplotypeModel(TreeModel treeModel,
			HaplotypeModel haplotypeModel) {
		
		double errorRate = 0;
		
		// clock model
		Parameter rateParameter = new Parameter.Default(StrictClockBranchRates.RATE, 1e-5, 0, 1);
		StrictClockBranchRates branchRateModel = new StrictClockBranchRates(rateParameter);
	
		Parameter freqs = new Parameter.Default("frequency", new double[]{0.25, 0.25, 0.25, 0.25});
		Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100.0);

		// Sub model
		FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
		HKY hky = new HKY(kappa, f);

		// siteModel
		GammaSiteModel siteModel = new GammaSiteModel(hky);
		Parameter mu = new Parameter.Default(
				GammaSiteModelParser.MUTATION_RATE, 1, 0, Double.POSITIVE_INFINITY);
		siteModel.setMutationRateParameter(mu);

		// Simulate halotypes
		if(errorRate>0){
			haplotypeModel.simulateSequence(errorRate, siteModel, hky, treeModel);
		}
		
		// treeLikelihood
		TreeLikelihoodExt treeLikelihood = new TreeLikelihoodExt(
				haplotypeModel, treeModel, siteModel, branchRateModel, null,
				false, false, true, false, false);
		treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);
		
		// treeLikelihood
//		TreeLikelihoodExt treeLikelihood = MCMCSetupHelper.setupTreeLikelihood(kappa, freqs,
//				haplotypeModel, treeModel, branchRateModel);
//	
		HashMap<String, Object> parameterList = new HashMap<String, Object>();
		parameterList.put("kappa", kappa);
		parameterList.put("freqs", freqs);
		parameterList.put("treeLikelihood", treeLikelihood);
		parameterList.put("branchRateModel", branchRateModel);
		
		return parameterList;
	}

	private static final double opTiny = 0.1;
	private static final double opSmall = 3;
	private static final double opMed = 10;
	private static final double opLarge = 30;
	private static final double opHuge = 100;
	
	public static TreeLikelihoodExt setupTreeLikelihood(Parameter kappa,
			Parameter freqs, HaplotypeModel haplotypeModel,
			TreeModel treeModel, StrictClockBranchRates branchRateModel) {
	
		// Sub model
		FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
		HKY hky = new HKY(kappa, f);
	
		// siteModel
		GammaSiteModel siteModel = new GammaSiteModel(hky);
		Parameter mu = new Parameter.Default(
				GammaSiteModelParser.MUTATION_RATE, 1, 0, Double.POSITIVE_INFINITY);
		siteModel.setMutationRateParameter(mu);
	
		// treeLikelihood
	
		TreeLikelihoodExt treeLikelihood = new TreeLikelihoodExt(
				haplotypeModel, treeModel, siteModel, branchRateModel, null,
				false, false, true, false, false);
		haplotypeModel.simulateSequence(treeLikelihood);
		treeLikelihood = new TreeLikelihoodExt(
				haplotypeModel, treeModel, siteModel, branchRateModel, null,
				false, false, true, false, false);
		
		treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);
		// SitePatternsExt patterns = new SitePatternsExt(haplotypeModel, null, 0, -1, 1, true);
		// TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, branchRateModel, null, false, false, true, false, false);
		
	
		return treeLikelihood;
	}

	public static List<MCMCOperator> defalutTreeOperators(OperatorSchedule schedule,
			TreeModel treeModel) {
	
			MCMCOperator operator;
//			List<MCMCOperator> OperatorList = new ArrayList<MCMCOperator>();
			
			
			Parameter allInternalHeights = treeModel.createNodeHeightsParameter(true, true, false);
			operator = new UpDownOperator(
					null,// new Scalable[] { new Scalable.Default(rateParameter) },
					new Scalable[] { new Scalable.Default(allInternalHeights) },
					0.75, 3.0, CoercionMode.COERCION_ON);
			operator.setWeight(opSmall);
			schedule.addOperator(operator);
	
			Parameter rootHeight = treeModel.getRootHeightParameter();
			rootHeight.setId("TREE_HEIGHT");
			operator = new ScaleOperator(rootHeight, 0.75);
			operator.setWeight(opSmall);
			schedule.addOperator(operator);
	
			Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);
			operator = new UniformOperator(internalHeights, 30.0);
			operator.setWeight(opMed);
			schedule.addOperator(operator);
	
			operator = new SubtreeSlideOperator(treeModel, 15.0, 1.0, true, false,
					false, false, CoercionMode.COERCION_ON);
			operator.setWeight(opMed);
			schedule.addOperator(operator);
	
			operator = new ExchangeOperator(ExchangeOperator.NARROW, treeModel, 15.0);
			operator.setWeight(opMed);
	//		schedule.addOperator(operator);
	
			operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 3.0);
			operator.setWeight(opSmall);
	//		schedule.addOperator(operator);
	
			operator = new WilsonBalding(treeModel, 3.0);
			operator.setWeight(opSmall);
			schedule.addOperator(operator);
	
			return null;
		}

	public static ArrayList<MCMCOperator> defalutOperators(OperatorSchedule schedule,
				HaplotypeModel haplotypeModel,
					Parameter... parameters) {
		
				MCMCOperator operator;
	//			ArrayList<MCMCOperator> OperatorList = new ArrayList<MCMCOperator>();
		
		//		operator = new SingleBaseOperator(haplotypeModel, 0);
		//		operator.setWeight(opSmall);
		//		schedule.addOperator(operator);
		//
		//		operator = new SingleBaseUniformOperator(haplotypeModel, 0);
		//		operator.setWeight(opSmall);
		//		schedule.addOperator(operator);
		
		//		operator = new SingleBaseEmpiricalOperator(haplotypeModel, 0);
		//		operator.setWeight(opLarge);
		//		schedule.addOperator(operator);
		
		//		operator = new SwapBasesMultiOperator(haplotypeModel, 12, CoercionMode.COERCION_ON);
		//		operator.setWeight(opLarge);
		//		schedule.addOperator(operator);
		//
		//		operator = new SwapBasesUniformOperator(haplotypeModel, 12, CoercionMode.COERCION_ON);
		//		operator.setWeight(opLarge);
		//		schedule.addOperator(operator);
		//
		//		operator = new SwapBasesEmpiricalOperator(haplotypeModel, 1, CoercionMode.COERCION_OFF);
		//		operator.setWeight(opLarge);
		//		schedule.addOperator(operator);
		
		//		operator = new SwapBasesEmpiricalOperator(haplotypeModel, 2, CoercionMode.COERCION_OFF);
		//		operator.setWeight(opLarge);
		//		schedule.addOperator(operator);
		//		operator = new SwapBasesEmpiricalOperator(haplotypeModel, 5, CoercionMode.COERCION_OFF);
		//		operator.setWeight(opLarge);
		//		schedule.addOperator(operator);
		//		operator = new SwapBasesEmpiricalOperator(haplotypeModel, 5, CoercionMode.COERCION_ON);
		//		operator.setWeight(opLarge);
		//		schedule.addOperator(operator);
		
//				operator = new HaplotypeRecombinationOperator(haplotypeModel, 12);
//				operator.setWeight(3.0);
		//		schedule.addOperator(operator);
		
				operator = new BaseSingleOperator(haplotypeModel);
				operator.setWeight(opLarge);
				schedule.addOperator(operator);
				
//				operator = new BasesMultiOperator(haplotypeModel, 3, CoercionMode.COERCION_ON);
//				operator.setWeight(opLarge/10.0);
//				schedule.addOperator(operator);
//				
//				operator = new HaplotypeSwapSectionOperator(haplotypeModel, 5, CoercionMode.COERCION_ON);
//				operator.setWeight(opLarge);
//				schedule.addOperator(operator);
//				
//				operator = new ColumnOperator(haplotypeModel, 7, CoercionMode.COERCION_OFF);
//				operator.setWeight(opLarge);
//				schedule.addOperator(operator);
				
				
				for (Parameter parameter : parameters) {
		
					String parameterName = parameter.getParameterName();
					
					if("kappa".equals(parameterName)){
						operator = new ScaleOperator(parameter, 0.75);
						operator.setWeight(opTiny);
						schedule.addOperator(operator);
					}
					else if("frequency".equals(parameterName)){
						operator = new DeltaExchangeOperator(parameter, new int[] { 1,
								1, 1, 1 }, 0.01, 0.1, false, CoercionMode.COERCION_ON);
						operator.setWeight(opTiny);
		//				operator.setWeight(opHuge);
						schedule.addOperator(operator);
						
//						operator = new BaseSingleFrequencyOperator(haplotypeModel, parameter);
//						operator.setWeight(opMed);
//						schedule.addOperator(operator);
//						

						
//						operator = new SwitchBaseFrequencyOperator(haplotypeModel, 0.8, 
//								parameter, CoercionMode.COERCION_ON);
//						operator.setWeight(opMed);
//						schedule.addOperator(operator);
						
						//good seq: low (switch) prob, most accepted with same posterior
						//			high switchProb, low accept, but with diff posterior
						//bad seq: high (switchProb), accept with different posterior
						//			low switch,  accept with same post
						
		//				operator = new SingleBaseFrequencyOperator(haplotypeModel, parameter);
		//				operator.setWeight(opMed);
		////				schedule.addOperator(operator);
					}
					else if("populationSize".equals(parameterName)){
						operator = new ScaleOperator(parameter, 0.75);
						operator.setWeight(opTiny);
						schedule.addOperator(operator);
					}
					
				}
				
				return null;
			}

	
}
