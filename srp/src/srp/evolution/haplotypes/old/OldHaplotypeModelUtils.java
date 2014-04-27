package srp.evolution.haplotypes.old;

import srp.evolution.haplotypes.AlignmentUtils;
import srp.evolution.shortreads.AlignmentMapping;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;

public class OldHaplotypeModelUtils {

	public static OldHaplotypeModel copyHaplotypeModel(
			OldHaplotypeModel haplotypeModel) {
		
		
        SimpleAlignment alignment = new SimpleAlignment();
        for (int j = 0; j < haplotypeModel.getHaplotypeCount(); j++) {
        	alignment.addSequence(haplotypeModel.getHaplotype(j));
		}

        OldHaplotypeModel copyHaplotypeModel = new OldHaplotypeModel(haplotypeModel.getAlignmentMapping(), alignment);
		return copyHaplotypeModel;
	}

	public static OldHaplotypeModel factory(Alignment shortReads, Alignment trueAlignment){
		
		AlignmentMapping alignmentMapping = new AlignmentMapping(shortReads);
		OldHaplotypeModel haplotypeModel = new OldHaplotypeModel(alignmentMapping, trueAlignment);
		return haplotypeModel;
	}

	public static OldHaplotypeModel createHaplotypeModel(String[] taxa_sequence){
		
		Alignment alignment = AlignmentUtils.createAlignment(taxa_sequence);
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(taxa_sequence);
		OldHaplotypeModel haplotypeModel = new OldHaplotypeModel(aMap, alignment);
		return haplotypeModel;
		
	}

	public static OldHaplotypeModel createHaplotypeModel(String[] shortRead, String[] taxa_sequence){
		
		Alignment alignment = AlignmentUtils.createAlignment(taxa_sequence);
		AlignmentMapping aMap = AlignmentUtils.createAlignmentMapping(shortRead);
		OldHaplotypeModel haplotypeModel = new OldHaplotypeModel(aMap, alignment);
		return haplotypeModel;
		
	}

	public static OldHaplotypeModel createHaplotypeModel(AlignmentMapping aMap, Alignment trueAlignment) {
		int noSeq = trueAlignment.getSequenceCount();
		OldHaplotypeModel haplotypeModel = new OldHaplotypeModel(aMap, noSeq);
		for (int i = 0; i < noSeq; i++) {
			String sequence = trueAlignment.getSequence(i).getSequenceString();
			haplotypeModel.getHaplotype(i).setSequenceString(sequence);
		}
		 
		return haplotypeModel;
	}

}
