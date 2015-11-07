package srp.operator.haplotypes;

import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.plaf.synth.SynthSeparatorUI;

import srp.dr.evolution.datatype.ShortReads;
import srp.evolution.OperationType;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class BaseDataSingleOperator extends AbstractSingleOperator {

	public static final String OPERATOR_NAME = BaseDataSingleOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.SINGLE;
	
	public BaseDataSingleOperator(HaplotypeModel haplotypeModel) {
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
		Haplotype haplotype = haplotypeModel.getHaplotype(hapIndex);
		
		ArrayList<int[]> listOfAvailableChar2 = haplotypeModel.getListOfAvailableChar2();
		

		int siteIndex;// = getNextSiteIndex();
		boolean isFix;// = false;
//		System.out.println(siteIndex +"\t"+ isFix);
		do {
			siteIndex = getNextSiteIndex();
//			isFix = srpMap.isFixSite(siteIndex);
			isFix = srpMap.isLowVarianceSite(siteIndex);
//			System.out.println("DO: "+siteIndex +"\t"+ isFix);
		} while (isFix);
//		System.out.println(siteIndex +"\t"+ isFix);
		
//		System.out.println();
		
//		char newChar = '-';
////		siteIndex = getNextSiteIndex();
//		boolean notFound = true;
//		do{
//			siteIndex = getNextSiteIndex();
////			int oldState = haplotype.getState(siteIndex);
//			int[] chars = listOfAvailableChar2.get(  siteIndex );
//			int size = chars.length;
//			
//			if(size > 1){
//				newChar = (char) chars[ MathUtils.nextInt(size) ];
////				if(newChar != DNA_CHARS[oldState]){
//					notFound = false;
////				}
//			}
////		else{
////			newChar = DNA_CHARS[ MathUtils.nextInt(4) ];
////		}
//		
//		
////		char newChar = getNextDiffBase(oldState);
////		char newChar = 'A';
//		} while(notFound);
		
//		newChar = getNextBase();
//		if(MathUtils.nextDouble()< 0.5){
//			newChar = cheat[hapIndex].charAt(siteIndex);
//		}
		
//		haplotypeModel.getGetMapToSrpAt(siteIndex);
//		newChar = (char) srpMap.nextBaseAt(siteIndex);
		char newChar = srpMap.nextBaseFreqAt(siteIndex);
		
		
//		System.out.println(newChar);
		haplotype.setCharAt(siteIndex, newChar);
		haplotypeModel.setOperationRecord(OP, hapIndex, siteIndex);

//		System.out.print( (ShortReads.NUCLEOTIDE_CHARS[oldState]) +"\t"+ newChar +"\ttrue"+ cheat[hapIndex].charAt(siteIndex)
//		+"\t"+ (cheat[hapIndex].charAt(siteIndex)==newChar) +"\t"
//		);
		
		haplotypeModel.endAlignmentModelOperation();

		return 0.0;
	}


	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}




}
