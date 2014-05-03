package srp.operator.haplotypes;

import srp.dr.evolution.datatype.ShortReads;
import srp.evolution.OperationType;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class BaseSingleOperator extends AbstractSingleOperator {

	public static final String OPERATOR_NAME = BaseSingleOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.SINGLE;
	
	public BaseSingleOperator(HaplotypeModel haplotypeModel) {
		super(haplotypeModel);
		
		site = 0;
		hap = 0;
	}

	static int site;
	static int hap;
	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startAlignmentModelOperation();

		int hapIndex = getNextHapIndex();
		int siteIndex = getNextSiteIndex();
		
//		if(hap== (haplotypeCount-1)){
//			hap =0;
////			System.out.println(hap +"a\ta"+ site);
//		}
//		if(site==(haplotypeLength-1)){
//			site = 0;
//			hap++;
////			System.out.println(hap +"\t"+ site);
//		}
//		site++;
////		System.out.println(hap +"\t"+ site);
//		hapIndex = hap;
//		siteIndex = site;

		Haplotype haplotype = haplotypeModel.getHaplotype(hapIndex);
		int oldState = haplotype.getState(siteIndex);
		char newChar = getNextDiffBase(oldState);
//		newChar = getNextBase();
//		if(MathUtils.nextDouble()< 0.5){
//			newChar = cheat[hapIndex].charAt(siteIndex);
//		}
		haplotype.setCharAt(siteIndex, newChar);
		haplotypeModel.setOperationRecord(OP, hapIndex, siteIndex);

//		System.out.print( (ShortReads.NUCLEOTIDE_CHARS[oldState]) +"\t"+ newChar +"\ttrue"+ cheat[hapIndex].charAt(siteIndex)
//		+"\t"+ (cheat[hapIndex].charAt(siteIndex)==newChar) +"\t"
//		);
		
		haplotypeModel.endAlignmentModelOperation();

		return 0.0;
	}
	static String cheat[] = new String[]{
		"AGTTAAGAGGCACAACTTTGGTGGATGTCAGTTGAGGTTGCTACTACACAAAAAGTACATACCTTACGGAACGTTCTCAATCACTGTGAGACGTTTGCCATGGTACACTCTGCGTCCACT",
		"AGTTAAGAGGCACAACTTTGGTGGATGTCAGTTGAGGTTGCTACTACACAAAAAGTACATACCTTACGGCACGTTCTCAATCACTGTGAGACGTTTGCCATGGTACACTCTGCGTCCACT",
		"AGTTAAGAGGCACAACTTTGGTGGATGTCAGTTGAGGTTGCTACTACACAAAAAGTACATACCTTACGGAACGTTCTCAATCACTGTGAGACGTTTGCCATGGTACACTCTGCGTCCACT",
		"AGCTAAGAGACACAACTATGGTGGATGTCAGTCAAGCCTGCAACTACATAGAAAGTGAATAAGTTACTGAATGCTCTCATTCACTGAGAGACGTTTACCATAATATACTAGGCGTCGACG",
		"AGCTAAGAGACACAACTATGGTGGATGTCAGTCAAGCCTGCAACTACATAGAAAGTGAATAAGCTACTGAATGCTCTCATTCACTGAGAGACGTTTACCATAATATACTAGGCGTCGACG",
		"AGCTAAGAGACACAACTATGGTGGATGTCAGTCAAGCCTGCAACTACATAGAAAGTGAATAAGTTACTGAATGCTCTCATTCACTGAGAGACGTTTACCATAATATACTTGGCGTCGACG",
		"AGCTAAGAGACACAACTATGGTGGATGTCAGTCAAGCCTCCTACTACATAGAAAGTGAATAAGTTACTGAATGCTCTCATTCACTGAGAGACGTTTGCCTTAATAGACTTGGCGTCGACG",
	};

	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}




}
