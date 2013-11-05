package srp.spectrum.operator;

import java.util.Arrays;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.Operation;
import srp.spectrum.SpectrumAlignmentModel;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.OperatorUtils;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

public abstract class AbstractSpectrumOperator extends AbstractCoercableOperator {

	public final static Operation OP = Operation.SWAPSINGLE;

	public SpectrumAlignmentModel spectrumModel;
	public int spectrumCount;
	public AlignmentMapping alignmentMapping;
	
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



