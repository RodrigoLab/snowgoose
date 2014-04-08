package srp.operator.haplotypes;

import srp.evolution.OperationRecord;
import srp.evolution.OperationType;
import srp.evolution.shortreads.AlignmentMapping;
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
	
	@Deprecated
	public AlignmentMapping alignmentMapping;

	public AbstractSingleOperator(HaplotypeModel haplotypeModel) {
		this.haplotypeModel = haplotypeModel;
		haplotypeCount = this.haplotypeModel.getHaplotypeCount();
		haplotypeLength = this.haplotypeModel.getHaplotypeLength();
		
		alignmentMapping = this.haplotypeModel.getAlignmentMapping();
	}

	

	public int getNextHapIndex(){
		return MathUtils.nextInt(haplotypeCount);
	}
	
	public int getNextSiteIndex(){
		return MathUtils.nextInt(haplotypeLength);
	}
	
	public char getNextBase(){
		int i = MathUtils.nextInt(DIMENSION);
		return DNA_CHARS[i];
	}



	@Override
	public String getPerformanceSuggestion() {
		return "";
	}
	
}