package srp.operator.haplotypes;

import java.util.Arrays;

import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.Operation;
import srp.shortreads.AlignmentMapping;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;

public abstract class AbstractBasesMultiOperator extends AbstractCoercableOperator {

	public final static Operation OP = Operation.SWAPMULTI;
	
	protected final int haplotypeLength;

	protected HaplotypeModel haplotypeModel;

	protected int swapLength;

	protected int[][] allPosChars;
	protected AlignmentMapping alignmentMapping;

	private double autoOptimize;
	private double scaleFactor;
	
	public AbstractBasesMultiOperator(HaplotypeModel haplotypeModel, int swapLength, CoercionMode mode) {
		super(mode);
		
		this.haplotypeModel = haplotypeModel;
		this.swapLength = swapLength;
		haplotypeLength = this.haplotypeModel.getHaplotypeLength();
		alignmentMapping = this.haplotypeModel.getAlignmentMapping();
		
		allPosChars = new int[2][haplotypeLength];
	
		
		checkParameterIsValid();
		
		scaleFactor = (int) (haplotypeLength*0.01);

		if (scaleFactor <1) {
			scaleFactor = 1;
		}
		
		convertToAutoOptimize(this.swapLength);
		
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
			swapLength =  1 + (int) Math.exp(autoOptimize*scaleFactor);
//			System.out.println(autoOptimize +"\t"+ Math.exp(autoOptimize*scaleFactor));
			
//			System.out.print("A=" + swapLength + "\t" + autoOptimize + "\t" +
//					"accept: " + getAcceptCount()/(double)getCount() + "\t"  );
			
			checkParameterIsValid();
			
//			System.out.print("newL:"+swapLength+" ");
	//		System.out.print("A\t" + swapFragmentLength + "\t" + autoOptimize + "\t"  );
	    }

	private double convertToAutoOptimize(int length) {
		swapLength = length;
		checkParameterIsValid();
		autoOptimize = Math.log(swapLength - 1)/scaleFactor;
	    return autoOptimize;
	}

	private void checkParameterIsValid() {
		if (swapLength > haplotypeLength){
			swapLength = haplotypeLength;
		}
	}
	
	protected void resetAllPosChars() {
		for (int i = 0; i < allPosChars.length; i++) {
			Arrays.fill(allPosChars[i], -1);
		}
	}


	@Override
	public double getRawParameter() {
		return swapLength;
	}
	

}
