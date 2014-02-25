package srp.spectrum.operator;

import java.util.HashSet;
import java.util.Set;

import com.google.common.primitives.Ints;

import srp.haplotypes.AlignmentMapping;
import srp.spectrum.SpectraParameter;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import dr.inference.model.Bounds;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.math.MathUtils;

public abstract class AbstractSpectrumOperator extends AbstractCoercableOperator {

	public static final int DIMENSION = SpectraParameter.DIMENSION;
    public static final Bounds<Double> BOUNDS = SpectraParameter.SPECTRA_BOUNDS;
    public static final double BOUNDS_LOWER = BOUNDS.getLowerLimit(0);
    public static final double BOUNDS_UPPER = BOUNDS.getUpperLimit(0);
    

    public SpectrumOperation OP;
    
    protected SpectrumAlignmentModel spectrumModel;
    @Deprecated
	protected AlignmentMapping alignmentMapping;
	
	protected final int spectrumCount;
	protected final int spectrumLength;


//	protected int swapLength;
//
//	protected int[][] allPosChars;
//
//	private double autoOptimize;
//	private double scaleFactor;
	
	public AbstractSpectrumOperator(SpectrumAlignmentModel spectrumModel, CoercionMode mode) {
		super(mode);
		this.spectrumModel = spectrumModel;
		spectrumCount = this.spectrumModel.getSpectrumCount();
		spectrumLength = this.spectrumModel.getSpectrumLength();
		alignmentMapping = this.spectrumModel.getAlignmentMapping();
	}


	public abstract SpectrumOperation getSpectrumOperation();

	
	public int[] generateSiteIndexs(int swapBasesCount, int spectrumLength) {

		Set<Integer> generated = new HashSet<Integer>();
		while (generated.size() < swapBasesCount)
		{
		    Integer next = MathUtils.nextInt(spectrumLength);
		    generated.add(next);
		}
		int[] siteIndexs = Ints.toArray(generated);

		return siteIndexs;
	}

//	public SpectrumOperation getSpectrumOperation() {
//		// TODO Auto-generated method stub
//		return OP;
//	}
	
//	public AbstractMultiBasesOperator(HaplotypeModel haplotypeModel, int swapLength, CoercionMode mode) {
//		super(mode);
//		
//		this.haplotypeModel = haplotypeModel;
//		this.swapLength = swapLength;
//		haplotypeLength = this.haplotypeModel.getHaplotypeLength();
//		alignmentMapping = this.haplotypeModel.getAlignmentMapping();
//		
//		allPosChars = new int[2][haplotypeLength];
//	
//		
//		checkParameterIsValid();
//		
//		scaleFactor = (int) (haplotypeLength*0.01);
//
//		if (scaleFactor <1) {
//			scaleFactor = 1;
//		}
//		
//		convertToAutoOptimize(this.swapLength);
//		
//	}

//	@Override
//	public double getCoercableParameter() {
//	    return autoOptimize;
//	}
//
//	@Override
//	public void setCoercableParameter(double autoOpt) {
//		convertFromAutoOptimizeToValue(autoOpt);
//	    
//	}

//	private void convertFromAutoOptimizeToValue(double autoOpt) {
//	    	autoOptimize = autoOpt;
//			swapLength =  1 + (int) Math.exp(autoOptimize*scaleFactor);
////			System.out.println(autoOptimize +"\t"+ Math.exp(autoOptimize*scaleFactor));
//			
////			System.out.print("A=" + swapLength + "\t" + autoOptimize + "\t" +
////					"accept: " + getAcceptCount()/(double)getCount() + "\t"  );
//			
//			checkParameterIsValid();
//			
////			System.out.print("newL:"+swapLength+" ");
//	//		System.out.print("A\t" + swapFragmentLength + "\t" + autoOptimize + "\t"  );
//	    }
//
//	private double convertToAutoOptimize(int length) {
//		swapLength = length;
//		checkParameterIsValid();
//		autoOptimize = Math.log(swapLength - 1)/scaleFactor;
//	    return autoOptimize;
//	}
//
//	private void checkParameterIsValid() {
//		if (swapLength > spectrumLength){
//			swapLength = spectrumLength;
//		}
//	}
//	
//	protected void resetAllPosChars() {
//		for (int i = 0; i < allPosChars.length; i++) {
//			Arrays.fill(allPosChars[i], -1);
//		}
//	}
//
//
//	@Override
//	public double getRawParameter() {
//		return swapLength;
//	}

//
//	@Override
//	public String getPerformanceSuggestion() {
//		return "";
//	}


}



