package srp.evolution.haplotypes;

import srp.dr.evolution.datatype.ShortReads;
import srp.evolution.shortreads.AlignmentMapping;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

public class AlignmentUtils {

	
	public static SimpleAlignment createAlignment(String[][] taxa_sequence){

        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(Nucleotides.INSTANCE);

        Taxon[] taxa = new Taxon[taxa_sequence[0].length]; // 6, 17
        System.out.println("Taxon len = " + taxa_sequence[0].length);
        System.out.println("Alignment len = " + taxa_sequence[1].length);
                         

        for (int i=0; i < taxa_sequence[0].length; i++) {
            taxa[i] = new Taxon(taxa_sequence[0][i].toString());

            Sequence sequence = new Sequence(taxa_sequence[1][i].toString());
            sequence.setTaxon(taxa[i]);
//            sequence.setDataType(dataType);

            alignment.addSequence(sequence);
        }
        return alignment;
    }
	

	public static SimpleAlignment createAlignment(String[] sequences){

        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(ShortReads.INSTANCE);

        Taxon[] taxa = new Taxon[sequences.length]; 
        for (int i=0; i < sequences.length; i++) {
            taxa[i] = new Taxon("taxa_"+i);

            Sequence sequence = new Sequence(sequences[i].toString());
            sequence.setTaxon(taxa[i]);
//            sequence.setDataType(dataType);

            alignment.addSequence(sequence);
        }
        
        return alignment;
        
    }
	
	public static AlignmentMapping createAlignmentMapping(String[] taxa_sequence){
		SimpleAlignment alignment = createAlignment(taxa_sequence);
		AlignmentMapping aMap = new AlignmentMapping(alignment);
		return aMap;
	}
}
