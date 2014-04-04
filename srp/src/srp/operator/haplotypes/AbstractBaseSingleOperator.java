package srp.operator.haplotypes;

import srp.evolution.haplotypes.old.OldHapOperation;
import srp.evolution.haplotypes.old.OldHaplotypeModel;
import srp.evolution.shortreads.AlignmentMapping;
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