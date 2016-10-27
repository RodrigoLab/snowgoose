package srp.operator.haplotypes;

import org.apache.commons.math3.util.FastMath;

import srp.evolution.OperationType;
import srp.evolution.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.math.MathUtils;

public abstract class AbstractMultiOperator extends AbstractHaplotypeOperator {

	public static final OperationType OP = OperationType.MULTI;

	private static final int MIN_BASE = 1;
	


	protected int basesCount;

	private double autoOptimize;
	private double scaleFactor;
	private int[] haplotypeLengthArray;
	
	public AbstractMultiOperator(HaplotypeModel haplotypeModel, int basesCount, CoercionMode mode) {
		super(haplotypeModel, mode);
		this.basesCount = basesCount;
		
//		allPosChars = new int[2][haplotypeLength];
		haplotypeLengthArray = new int[haplotypeLength];
		for (int i = 0; i < haplotypeLengthArray.length; i++) {
			haplotypeLengthArray[i] = i;
		}
		
		
		checkParameterIsValid();
		
		scaleFactor = (int) (haplotypeLength*0.01);

		if (scaleFactor <1) {
			scaleFactor = 1;
		}
		
		convertToAutoOptimize(this.basesCount);
		
	}


	
//	public int getNextSiteIndex(){
//		return MathUtils.nextInt(haplotypeLength);
//	}
	
	public int[] generateUniqueSites(int m) {
		int[] siteIndexs = randomSampleSites(m, haplotypeLengthArray);
		return siteIndexs;
	}

	
	public static int[] randomSampleSites(int m, int[] sampleArray){ //time:0.5
		int[] sites = new int[m];
		int length = sampleArray.length;
		for (int i = 0; i < m; i++) {
	        int pos = i + MathUtils.nextInt(length - i);
//	        T tmp = items.get(pos);
//	        items.set(pos, items.get(i));
//	        items.set(i, tmp);
	        int tmp = sampleArray[pos];
	        sampleArray[pos] = sampleArray[i];
	        sampleArray[i] = tmp;
	        sites[i] = tmp;
	        
	    }
//	    return items.subList(0, m);
	    return sites;
	}
	
	public OperationType getOperationType() {
		return OP;
	}



	@Override
	public final String getPerformanceSuggestion() {
		String s = "Tuning "+basesCount; 
		return s;
	
	}



	@Override
	public double getCoercableParameter() {
	    return autoOptimize;
	}

	@Override
	public void setCoercableParameter(double autoOpt) {
		convertFromAutoOptimizeToValue(autoOpt);
	}

	@Override
	public double getRawParameter() {
        return basesCount;
    }

    private void convertFromAutoOptimizeToValue(double autoOpt) {
		autoOptimize = autoOpt;
		basesCount = MIN_BASE + (int) FastMath.exp(autoOptimize);
	
		checkParameterIsValid();
	}



	private double convertToAutoOptimize(int length) {
		basesCount = length;
		checkParameterIsValid();
		autoOptimize = Math.log(basesCount - MIN_BASE);
	    return autoOptimize;
	}



	private void checkParameterIsValid() {
		if (basesCount >= haplotypeLength){
			basesCount = haplotypeLength-1;
		}
	}



	@Override
	public double getTargetAcceptanceProbability() {
	    return 0.234;
	}



	@Override
	public String toString() {
	    return getOperatorName() + "(AbstractToString())";
	}

	
}
