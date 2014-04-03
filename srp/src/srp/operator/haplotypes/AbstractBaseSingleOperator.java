package srp.operator.haplotypes;

import srp.haplotypes.old.OldHaplotypeModel;
import srp.haplotypes.old.OldHapOperation;
import srp.shortreads.AlignmentMapping;
import dr.inference.operators.SimpleMCMCOperator;

public abstract class AbstractBaseSingleOperator extends SimpleMCMCOperator {

	public final static OldHapOperation OP = OldHapOperation.SWAPSINGLE;

	public OldHaplotypeModel haplotypeModel;
	public int haplotypeCount;
	public AlignmentMapping alignmentMapping;

	public AbstractBaseSingleOperator(OldHaplotypeModel haplotypeModel) {
		this.haplotypeModel = haplotypeModel;
		haplotypeCount = this.haplotypeModel.getHaplotypeCount();
		alignmentMapping = this.haplotypeModel.getAlignmentMapping();
	}

	@Override
	public String getPerformanceSuggestion() {
		return "";
	}


}