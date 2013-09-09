package srp.haplotypes.operator;

import java.util.Arrays;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.Operation;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;

public abstract class AbstractSwapBasesOperator  extends AbstractCoercableOperator {

	public final static Operation OP = Operation.SWAPMULTI;
	
	protected final int haplotypeLength;
	protected final int haplotypeCount;

	protected HaplotypeModel haplotypeModel;
	
	protected int swapLength;
	protected double autoOptimize;
	protected double scaleFactor;

	protected int[][] allPosChars;
	protected AlignmentMapping alignmentMapping; 
	
	public AbstractSwapBasesOperator(HaplotypeModel haplotypeModel, int swapLength, CoercionMode mode) {
		super(mode);
		
		this.haplotypeModel = haplotypeModel;
		haplotypeLength = this.haplotypeModel.getHaplotypeLength();
		haplotypeCount  = this.haplotypeModel.getHaplotypeCount();
		alignmentMapping = this.haplotypeModel.getAlignmentMapping();
		
		allPosChars = new int[2][haplotypeLength];
		
		if (swapLength < 1){
    		swapLength = 1;
    	}
		this.swapLength = swapLength;

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
	

}
