package srp.haplotypes.operator;

import srp.haplotypes.HaplotypeModel;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;

public abstract class AbstractHaplotypeOperator  extends AbstractCoercableOperator {

	protected final int haplotypeLength;
	protected final int haplotypeCount;

	protected HaplotypeModel haplotypeModel;
	
	protected int swapLength;
	protected double autoOptimize;
	protected double scaleFactor;

	
	
	public AbstractHaplotypeOperator(HaplotypeModel haplotypeModel, int swapLength, CoercionMode mode) {
		super(mode);
		
		this.haplotypeModel = haplotypeModel;
		this.haplotypeLength = this.haplotypeModel.getHaplotypeLength();
		this.haplotypeCount  = this.haplotypeModel.getHaplotypeCount();
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

}
