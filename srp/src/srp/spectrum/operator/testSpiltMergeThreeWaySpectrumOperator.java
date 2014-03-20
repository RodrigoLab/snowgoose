package srp.spectrum.operator;

import java.util.Arrays;

import javax.swing.text.TabableView;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.util.FastMath;

import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class testSpiltMergeThreeWaySpectrumOperator extends AbstractSpectrumOperator {

	public static final String OPERATOR_NAME = testSpiltMergeThreeWaySpectrumOperator.class.getSimpleName();
	public static final SpectrumOperation OP = SpectrumOperation.DELTA_MULTI;
	
    private final int[] parameterWeights;
    private double delta;
    private int swapBasesCount;
    private double autoOptimize;
    
//    private double[] debugList = new double[8];
//	private int scaleFactor=1;
//	public static double totalCount = 0;
	
	public testSpiltMergeThreeWaySpectrumOperator(SpectrumAlignmentModel spectrumModel, 
			double delta, int baseCount, CoercionMode mode) {
		super(spectrumModel, mode);
		
		
		this.delta = delta;
		this.swapBasesCount = baseCount;
        setWeight(1.0);

        parameterWeights = new int[this.spectrumModel.getDataType().getStateCount()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }
        convertToAutoOptimize(this.swapBasesCount);

	}

	public static final double MAX = 0.97;
	public static final double HALF = 0.49;
	public static final double THIRD = 0.33;
	public static final double MIN = 0.01;
	@Override
	public double doOperation() throws OperatorFailedException {

		spectrumModel.startSpectrumOperation();

//		spectrumModel.swapHaplotypeSingleBase(OP, posChar);
		int spectrumIndex = MathUtils.nextInt(spectrumCount);
        SpectraParameter[] spectra = new SpectraParameter[swapBasesCount];

        Spectrum spectrum = spectrumModel.getSpectrum(spectrumIndex);
		int[] siteIndexs = generateUniqueSites(swapBasesCount);
		

		for (int i = 0; i < swapBasesCount; i++) {
			
			spectra[i] = spectrum.getSpectra(siteIndexs[i]);
			
			int maxIndex = -1;
			int[] equalIndexs = new int[2];
			int[] equalTriple = new int[3];
			int count = 0;
			double[] oldFreq = new double[DIMENSION];
			for (int j = 0; j < DIMENSION; j++) {
				oldFreq[j] = spectra[i].getFrequency(j);
				if(oldFreq[j] == MAX ){
					maxIndex = j;
				}
				if(oldFreq[j] == HALF ){
					equalIndexs[count++] = j;
				}
				if(oldFreq[j] == THIRD ){
					equalTriple[count++] = j;
				}
			}
			
			if(maxIndex!= -1){//split
				int dim2 = getAnotherDimension(maxIndex);
				boolean split = MathUtils.nextBoolean();
				if(split){
					spectra[i].setFrequency(maxIndex, HALF);
					spectra[i].setFrequency(dim2, HALF);
				}
				else{
					int dim3;
					do {
						dim3 = MathUtils.nextInt(DIMENSION);
					} while (dim3 == maxIndex || dim3 == dim2);
					spectra[i].setFrequency(maxIndex, THIRD);
					spectra[i].setFrequency(dim2, THIRD);
					spectra[i].setFrequency(dim3, THIRD);
				}
//				System.out.println("S "+Arrays.toString(spectra[i].getFrequencies()));
			}
			
			else{//merge
//				System.out.println(Arrays.toString(spectra[i].getFrequencies()));
				if(count==2){//2 to 1
					boolean mergeOrder = MathUtils.nextBoolean();
					if(mergeOrder){
						spectra[i].setFrequency(equalIndexs[0], MIN);
						spectra[i].setFrequency(equalIndexs[1], MAX);
					}
					else{
						spectra[i].setFrequency(equalIndexs[0], MAX);
						spectra[i].setFrequency(equalIndexs[1], MIN);
					}
				}
				else if(count==3){//3 to 1
					int max = MathUtils.nextInt(count);
					for (int j = 0; j < count; j++) {
						if(j==max){
							spectra[i].setFrequency(equalTriple[j], MAX);
						}
						else{
							spectra[i].setFrequency(equalTriple[j], MIN);
						}
					}
//					System.out.print(max +"\t" );
				}
//				System.out.println(Arrays.toString(spectra[i].getFrequencies())+"\n");
			}
		}

		spectrumModel.setSpectrumOperationRecord(OP, spectrumIndex, siteIndexs);
		
		spectrumModel.endSpectrumOperation();
		
		return 0;
	}


	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}


//    @Override
//	public double getCoercableParameter() {
//        return Math.log(delta);
//    }
//
//    @Override
//	public void setCoercableParameter(double value) {
//        delta = Math.exp(value);
//
//    }

	@Override
	public double getCoercableParameter() {
	    return autoOptimize;
	}

	@Override
	public void setCoercableParameter(double autoOpt) {
		convertFromAutoOptimizeToValue(autoOpt);
	}

	private void convertFromAutoOptimizeToValue(double autoOpt) {
	    	autoOptimize = autoOpt;
	    	swapBasesCount =  MIN_BASE + (int) FastMath.exp(autoOptimize);
//			System.out.println(autoOptimize +"\t"+ Math.exp(autoOptimize*scaleFactor));
			
//			System.out.print("A=" + swapLength + "\t" + autoOptimize + "\t" +
//					"accept: " + getAcceptCount()/(double)getCount() + "\t"  );
			
			checkParameterIsValid();
			
//			System.out.print("newL:"+swapLength+" ");
	//		System.out.print("A\t" + swapFragmentLength + "\t" + autoOptimize + "\t"  );
	    }

	private double convertToAutoOptimize(int length) {
		swapBasesCount = length;
		checkParameterIsValid();
		autoOptimize = Math.log(swapBasesCount - MIN_BASE);
	    return autoOptimize;
	}

	private void checkParameterIsValid() {
		if (swapBasesCount > spectrumLength){
			swapBasesCount = spectrumLength;
		}
	}
	
    @Override
	public double getRawParameter() {
//        return delta;
        return swapBasesCount;
    }

    @Override
	public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    @Override
	public final String getPerformanceSuggestion() {
    	String s = "Tuning "+delta; 
    	return s;

    }

    @Override
	public String toString() {
        return getOperatorName() + "(windowsize=" + delta + ")";
    }


	@Override
	public SpectrumOperation getSpectrumOperation() {
		return OP;
	}

	
}





