package srp.operator.haplotypes;

import srp.evolution.OperationType;
import srp.haplotypes.HaplotypeModel;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

public abstract class AbstractSingleOperator extends SimpleMCMCOperator {

	public static final DataType DATATYPE = Nucleotides.INSTANCE;
	public static final OperationType OP = OperationType.SINGLE;
	public static final char[] DNA_CHARS = {'A','C','G','T'};
	public static final int DIMENSION = DNA_CHARS.length;
	
	public final int haplotypeCount;
	public final int haplotypeLength;
	 
	protected HaplotypeModel haplotypeModel;
	
	public AbstractSingleOperator(HaplotypeModel haplotypeModel) {
		this.haplotypeModel = haplotypeModel;
		haplotypeCount = this.haplotypeModel.getHaplotypeCount();
		haplotypeLength = this.haplotypeModel.getHaplotypeLength();
		
	}

	

	public int getNextHapIndex(){
		return MathUtils.nextInt(haplotypeCount);
	}
	
	public int getNextSiteIndex(){
		return MathUtils.nextInt(haplotypeLength);
	}
	
	public char getNextDiffBase(int oldState){
		int i = oldState;
		do {
			i = MathUtils.nextInt(DIMENSION);
		} while (i == oldState);

		return DNA_CHARS[i];
	}



	@Override
	public String getPerformanceSuggestion() {
		return "";
	}
	
}