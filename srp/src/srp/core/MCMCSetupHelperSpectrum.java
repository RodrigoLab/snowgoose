package srp.core;

import java.util.HashMap;

import org.apache.commons.math3.random.RandomDataGenerator;

import srp.haplotypes.HaplotypeModel;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.operator.DeltaExchangeColumnSpectrumOperator;
import srp.spectrum.operator.DeltaExchangeMultiSpectrumOperator;
import srp.spectrum.operator.DeltaExchangeSingleSpectrumOperator;
import srp.spectrum.operator.DirichletAlphaSpectrumOperator;
import srp.spectrum.operator.DirichletSpectrumOperator;
import srp.spectrum.operator.RecombinationSpectrumOperator;
import srp.spectrum.operator.RecombineSectionSpectrumOperator;
import srp.spectrum.operator.SwapColumnSpectrumOperator;
import srp.spectrum.operator.SwapMultiSpectrumOperator;
import srp.spectrum.operator.SwapSingleSpectrumOperator;
import srp.spectrum.operator.SwapSubColumnSpectrumOperator;
import srp.spectrum.treelikelihood.SpectrumTreeLikelihood;
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
import dr.ext.TreeLikelihoodExt;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.DeltaExchangeOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.Scalable;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.UniformOperator;
import dr.inference.operators.UpDownOperator;

public class MCMCSetupHelperSpectrum extends MCMCSetupHelper {

	private static final double opTiny = 0.1;
	private static final double opSmall = 3;
	private static final double opMed = 15;
	private static final double opLarge = 20;
	private static final double opSpectrum =150;
	
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
			parameterList.put("rate", rateParameter);
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
//		schedule.addOperator(operator);

		Parameter rootHeight = treeModel.getRootHeightParameter();
		rootHeight.setId("TREE_HEIGHT");
		operator = new ScaleOperator(rootHeight, 0.75);
		operator.setWeight(opSmall);
		schedule.addOperator(operator);

		Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);
		operator = new UniformOperator(internalHeights, 30.0);
		operator.setWeight(opLarge);
		schedule.addOperator(operator);

		operator = new SubtreeSlideOperator(treeModel, 15.0, 1.0, true, false,
				false, false, CoercionMode.COERCION_ON);
		operator.setWeight(opMed);
//		schedule.addOperator(operator);

		operator = new ExchangeOperator(ExchangeOperator.NARROW, treeModel, 15.0);
		operator.setWeight(opMed);
//		schedule.addOperator(operator);

		operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 3.0);
		operator.setWeight(opSmall);
//		schedule.addOperator(operator);

		operator = new WilsonBalding(treeModel, 3.0);
		operator.setWeight(opSmall);
//		schedule.addOperator(operator);

		return schedule;
	}

	public static OperatorSchedule defalutSpectrumOperators(
			OperatorSchedule schedule,
			SpectrumAlignmentModel spectrumModel, Parameter... parameters) {

//COERCION_OFF
//		Operator                                          Tuning   Count      Time     Time/Op  Pr(accept) 
//0.2		SingleSpectrumDeltaExchangeOperator                       667156     75941    0.11     0.0565      
//0.1		ColumnSpectrumDeltaExchangeOperator                       3332844    234268   0.07     0.0209   
		
		
		MCMCOperator operator;
		
		operator = new DeltaExchangeSingleSpectrumOperator(spectrumModel, 0.2, CoercionMode.COERCION_OFF);
//		0.05 delta, accept 0.8
		operator.setWeight(opSpectrum);
//		schedule.addOperator(operator);
//		Operator                                          Tuning   Count      Time     Time/Op  Pr(accept) 
//		SingleSpectrumDeltaExchangeOperator               0.205   1000000    117482   0.12     0.2375      



		
		operator = new DeltaExchangeColumnSpectrumOperator(spectrumModel, 0.05, CoercionMode.COERCION_OFF);
		operator.setWeight(opSpectrum*10);//fast! 
//		schedule.addOperator(operator);
//		COERCION_ON
//		Operator                                          Tuning   Count      Time     Time/Op  Pr(accept) 
//		ColumnSpectrumDeltaExchangeOperator               0.012   1000000    141392   0.14     0.2344 
		
		
		operator = new DeltaExchangeMultiSpectrumOperator(spectrumModel, 0.1,
				6, CoercionMode.COERCION_OFF);
		operator.setWeight(opSpectrum*2);//fix delta at 0.02(~7bp) or 0.05(3bp)
//		schedule.addOperator(operator);
		
		operator = new RecombinationSpectrumOperator(spectrumModel);
		operator.setWeight(opSpectrum/10); //(3bp)
//		schedule.addOperator(operator);


		operator = new RecombineSectionSpectrumOperator(spectrumModel, 6,
				CoercionMode.COERCION_ON);
		operator.setWeight(opMed); //(3bp)
		schedule.addOperator(operator);

//		Operator analysis
//		Operator                                          Tuning   Count      Time     Time/Op  Pr(accept)  Performance suggestion
//		SingleSpectrumDeltaExchangeOperator                       67009      8749     0.13     0.0076      low	Tuning 0.2
//		ColumnSpectrumDeltaExchangeOperator                       666179     14209    0.02     0.0011      very low	Tuning 0.05
//		RecombinationSpectrumOperator                             133587     1464373  10.96    0.0425      low	Tuning 12
//		SwapSectionSpectrumOperator                               133225     63241    0.47     0.0419      low	Tuning 12

		operator = new SwapSingleSpectrumOperator(spectrumModel, false);
		operator.setWeight(opSpectrum*1);
//		schedule.addOperator(operator);
		
		operator = new SwapMultiSpectrumOperator(spectrumModel, 3, CoercionMode.COERCION_ON);
		operator.setWeight(opSpectrum*1);
//		schedule.addOperator(operator);
		
		operator = new SwapColumnSpectrumOperator(spectrumModel, CoercionMode.COERCION_OFF);
		operator.setWeight(opSpectrum*5);
//		schedule.addOperator(operator);
		
		operator = new SwapSubColumnSpectrumOperator(spectrumModel, 1, CoercionMode.COERCION_OFF);
		operator.setWeight(opSpectrum*1);
//		schedule.addOperator(operator);
		
		operator = new SwapSubColumnSpectrumOperator(spectrumModel, 3, CoercionMode.COERCION_ON);
		operator.setWeight(opSpectrum*1);
//		schedule.addOperator(operator);

		// SwapSubColumnSpectrumOperator 39971 12299 0.31 0.2589 good Tuning 1
		// SwapSubColumnSpectrumOperator 40076 13092 0.33 0.0485 low Tuning 2
		// SwapSubColumnSpectrumOperator 39805 13713 0.34 0.0071 low Tuning 3
		// SwapSubColumnSpectrumOperator 40084 14282 0.36 0.0 very low Tuning 4
		// SwapSubColumnSpectrumOperator 2.0 40054 12469 0.31 0.2342 good Tuning 2

		operator = new DirichletSpectrumOperator(spectrumModel, 3, CoercionMode.COERCION_ON);
		operator.setWeight(opSpectrum*1);
		schedule.addOperator(operator);
		operator = new DirichletSpectrumOperator(spectrumModel, 1, CoercionMode.COERCION_OFF);
		operator.setWeight(opSpectrum*1);
		schedule.addOperator(operator);
//		operator = new DirichletSpectrumOperator(spectrumModel, 12, CoercionMode.COERCION_OFF);
//		operator.setWeight(opSpectrum*1);
//		schedule.addOperator(operator);

		operator = new DirichletAlphaSpectrumOperator(spectrumModel, 100, CoercionMode.COERCION_ON);
		operator.setWeight(opSpectrum*1);
		schedule.addOperator(operator);
//		operator = new DirichletAlphaSpectrumOperator(spectrumModel, 10, CoercionMode.COERCION_OFF);
//		operator.setWeight(opSpectrum*1);
//		schedule.addOperator(operator);
//		operator = new DirichletAlphaSpectrumOperator(spectrumModel, 1000, CoercionMode.COERCION_OFF);
//		operator.setWeight(opSpectrum*1);
//		schedule.addOperator(operator);

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
				schedule.addOperator(operator);
				
			}
			else if("populationSize".equals(parameterName)){
				operator = new ScaleOperator(parameter, 0.75);
				operator.setWeight(opSmall);
				schedule.addOperator(operator);
			}
			
		}
		
		return schedule;

	}
	
}
