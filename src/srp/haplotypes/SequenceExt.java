package srp.haplotypes;

import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

public class SequenceExt extends Sequence {
	
	StringBuilder sequenceString = new StringBuilder();
	public SequenceExt() {
		// TODO Auto-generated constructor stub

	}

	public SequenceExt(String sequence) {
		super(sequence);
		// TODO Auto-generated constructor stub
	}

	public SequenceExt(Sequence sequence) {
		super(sequence);
		// TODO Auto-generated constructor stub
	}

	public SequenceExt(Taxon taxon, String sequence) {
		super(taxon, sequence);
		// TODO Auto-generated constructor stub
	}

}
