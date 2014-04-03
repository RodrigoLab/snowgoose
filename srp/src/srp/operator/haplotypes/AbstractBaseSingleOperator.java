package srp.operator.haplotypes;

import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.Operation;
import srp.shortreads.AlignmentMapping;
import dr.inference.operators.SimpleMCMCOperator;

public abstract class AbstractBaseSingleOperator extends SimpleMCMCOperator {

	public final static Operation OP = Operation.SWAPSINGLE;

	public HaplotypeModel haplotypeModel;
	public int haplotypeCount;
	public AlignmentMapping alignmentMapping;

	public AbstractBaseSingleOperator(HaplotypeModel haplotypeModel) {
		this.haplotypeModel = haplotypeModel;
		haplotypeCount = this.haplotypeModel.getHaplotypeCount();
		alignmentMapping = this.haplotypeModel.getAlignmentMapping();
	}

	@Override
	public String getPerformanceSuggestion() {
		return "";
	}


}