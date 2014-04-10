package srp.operator.haplotypes;

import org.apache.commons.math3.util.FastMath;

import srp.evolution.OperationType;
import srp.haplotypes.HaplotypeModel;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.math.MathUtils;

public abstract class AbstractMultiOperator extends AbstractHaplotypeOperator {

	public static final OperationType OP = OperationType.MULTI;
	public static final DataType DATATYPE = Nucleotides.INSTANCE;
	public static final char[] DNA_CHARS = {'A','C','G','T'};
	public static final int DIMENSION = DNA_CHARS.length;
	
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

	public int getNextHapIndex(){
		return MathUtils.nextInt(haplotypeCount);
	}
	
//	public int getNextSiteIndex(){
//		return MathUtils.nextInt(haplotypeLength);
//	}
	
	public char getNextBase(){
		int i = MathUtils.nextInt(DIMENSION);
		return DNA_CHARS[i];
	}
	
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
		if (basesCount > haplotypeLength){
			basesCount = haplotypeLength;
		}
	}
	
    @Override
	public double getRawParameter() {
        return basesCount;
    }

    @Override
	public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    @Override
	public final String getPerformanceSuggestion() {
    	String s = "Tuning "+basesCount; 
    	return s;

    }

    @Override
	public String toString() {
        return getOperatorName() + "(AbstractToString())";
    }


	
	public OperationType getOperationType() {
		return OP;
	}

	
}
