package srp.core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import javax.swing.JPopupMenu.Separator;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.MathUtils;

import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.operator.ColumnOperator;
import srp.haplotypes.operator.HaplotypeRecombinationOperator;
import srp.haplotypes.operator.HaplotypeSwapSectionOperator;
import srp.haplotypes.operator.SingleBaseEmpiricalOperator;
import srp.haplotypes.operator.SingleBaseFrequencyOperator;
import srp.haplotypes.operator.SingleBaseOperator;
import srp.haplotypes.operator.SingleBaseUniformOperator;
import srp.haplotypes.operator.MultiBasesEmpiricalOperator;
import srp.haplotypes.operator.MultiBasesOperator;
import srp.haplotypes.operator.MultiBasesUniformOperator;
import srp.haplotypes.operator.SwitchBaseFrequencyOperator;
import srp.rj.operator.RJTreeOperator;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.likelihood.SpectrumTreeLikelihood;
import srp.spectrum.operator.ColumnSpectrumDeltaExchangeOperator;
import srp.spectrum.operator.MultiSpectrumDeltaExchangeOperator;
import srp.spectrum.operator.RecombinationSpectrumOperator;
import srp.spectrum.operator.SingleSpectrumDeltaExchangeOperator;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.ext.TreeLikelihoodExt;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.DeltaExchangeOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.Scalable;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.operators.UniformOperator;
import dr.inference.operators.UpDownOperator;
import dr.rj.RJTreeModel;

public class MCMCSetupHelperSpectrum extends MCMCSetupHelper {

	private static final double opTiny = 0.1;
	private static final double opSmall = 3;
	private static final double opMed = 12;
	private static final double opLarge = 50;
	private static final double opSpectrum = 24;
	
	public static HashMap<String, Object> setupSpectrumTreeLikelihoodSpectrumModel(
				TreeModel treeModel, SpectrumAlignmentModel spectrumModel) {
		RandomDataGenerator r = new RandomDataGenerator();
//		r.nextSample(c, k)
	//		double errorRate = 0.01;
			
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
	
	//		// Simulate halotypes
	//		if(errorRate>0){
	//			haplotypeModel.simulateSequence(errorRate, siteModel, hky, treeModel);
	//		}
			
			// treeLikelihood
			SpectrumTreeLikelihood treeLikelihood = new SpectrumTreeLikelihood(
					spectrumModel, treeModel, siteModel, branchRateModel, 
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

	public static OperatorSchedule defalutTreeOperators(
			OperatorSchedule schedule, TreeModel treeModel) {

		MCMCOperator operator;
//		List<MCMCOperator> OperatorList = new ArrayList<MCMCOperator>();
		
		
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
		schedule.addOperator(operator);

		operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 3.0);
		operator.setWeight(opSmall);
		schedule.addOperator(operator);

		operator = new WilsonBalding(treeModel, 3.0);
		operator.setWeight(opSmall);
		schedule.addOperator(operator);

		return schedule;
	}

	public static OperatorSchedule defalutSpectrumOperators(
			OperatorSchedule schedule,
			SpectrumAlignmentModel spectrumModel, Parameter... parameters) {


		MCMCOperator operator;
		
		operator = new SingleSpectrumDeltaExchangeOperator(spectrumModel, 0.2, CoercionMode.COERCION_ON);
//		0.05 delta, accept 0.8
		operator.setWeight(opSpectrum);
//		schedule.addOperator(operator);
		
		
		operator = new ColumnSpectrumDeltaExchangeOperator(spectrumModel, 0.05, CoercionMode.COERCION_OFF);
		operator.setWeight(opSpectrum);
//		schedule.addOperator(operator);
		
		operator = new MultiSpectrumDeltaExchangeOperator(spectrumModel, 0.05, CoercionMode.COERCION_OFF);
		operator.setWeight(opSpectrum);
		schedule.addOperator(operator);
		
		operator = new RecombinationSpectrumOperator(spectrumModel, 10, CoercionMode.COERCION_OFF);
		operator.setWeight(opSpectrum);
		schedule.addOperator(operator);
		
		for (Parameter parameter : parameters) {

			String parameterName = parameter.getParameterName();
			
			if("kappa".equals(parameterName)){
				operator = new ScaleOperator(parameter, 0.75);
				operator.setWeight(opTiny);
//				schedule.addOperator(operator);
			}
			else if("frequency".equals(parameterName)){
				operator = new DeltaExchangeOperator(parameter, new int[] { 1,
						1, 1, 1 }, 0.01, 0.1, false, CoercionMode.COERCION_ON);
				operator.setWeight(opTiny);
//				schedule.addOperator(operator);
				
			}
			else if("populationSize".equals(parameterName)){
				operator = new ScaleOperator(parameter, 0.75);
				operator.setWeight(opSmall);
//				schedule.addOperator(operator);
			}
			
		}
		
		return schedule;

	}

	
}
