package srp.operator.haplotypes;

import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.plaf.synth.SynthSeparatorUI;

import org.omg.CORBA.DATA_CONVERSION;

import srp.dr.evolution.datatype.ShortReads;
import srp.evolution.OperationType;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import srp.evolution.shortreads.SiteType;
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
		
// isLowVarianceSite
		
		
		int siteIndex;// = getNextSiteIndex();
		boolean isFix;// = false;
//		System.out.println(siteIndex +"\t"+ isFix);
		char newChar = 0;
		double logq = 0;
		
//		if(MathUtils.nextDouble()< 0){
//			do {
//				siteIndex = getNextSiteIndex();
//	//			isFix = srpMap.isFixSite(siteIndex);
//				isFix = srpMap.isLowVarianceSite(siteIndex);
//	//			System.out.println("DO: "+siteIndex +"\t"+ isFix);
//			} while (isFix);
////			System.out.println(siteIndex +"\t"+ isFix);
//			newChar = srpMap.nextBaseFreqAt(siteIndex);
//		}
//		else{
//			do {
//				siteIndex = getNextSiteIndex();
//	//			isFix = srpMap.isFixSite(siteIndex);
//				isFix = srpMap.isFixSite(siteIndex);
//	//			System.out.println("DO: "+siteIndex +"\t"+ isFix);
//			} while (!isFix);
//			
//			int oldState = haplotype.getState(siteIndex);
//			newChar = getNextDiffBase(oldState);
//		}
//		siteIndex = getNextSiteIndex();
//		int oldState = haplotype.getState(siteIndex);
//		newChar = getNextBase();
//		newChar = getNextDiffBase(oldState);
//		newChar = srpMap.nextBaseFreqAt(siteIndex);
		
		SiteType siteType;
//		do {
//			siteIndex = getNextSiteIndex();
////			isFix = srpMap.isFixSite(siteIndex);
//			siteType = srpMap.getSiteType(siteIndex);
////			System.out.println("DO: "+siteIndex +"\t"+ isFix);
//		} while (siteType == SiteType.FIXED);

		
//		switch (siteType) {
//		case LOW_VAR:
//			int oldState = haplotype.getState(siteIndex);
////			newChar = getNextDiffBase(oldState);		
//			newChar = srpMap.nextBaseFreqAt(siteIndex);
//			break;
//		case HIGH_VAR:
//			newChar = srpMap.nextBaseFreqAt(siteIndex);
//			break;
//		default:
//			System.exit(-1);
//			break;
//		}
	
		if(MathUtils.nextDouble()< 0.50){
			do {
				siteIndex = getNextSiteIndex();
				siteType = srpMap.getSiteType(siteIndex);
				
				newChar = srpMap.nextBaseFreqAt(siteIndex);
				char oldChar = haplotype.getChar(siteIndex);
				double oldCharProb = srpMap.getFreqAtSiteChar(siteIndex, oldChar);
				double newCHarProb = srpMap.getFreqAtSiteChar(siteIndex, newChar);
				logq = Math.log(newCHarProb/oldCharProb);
//				System.out.println(siteType +"\t"+ oldChar +"\t"+ newChar +"\t"+ logq +"\t"+ 
//						Arrays.toString(srpMap.getFreqAtSite(siteIndex))	);
//				int oldState = haplotype.getState(siteIndex);
//				newChar = getNextDiffBase(oldState);
			} while (siteType != SiteType.HIGH_VAR);
		}
		
		else{
//			do {
				siteIndex = getNextSiteIndex();
				siteType = srpMap.getSiteType(siteIndex);
				int oldState = haplotype.getState(siteIndex);
				newChar = getNextDiffBase(oldState);
				
//				double[] charProb = srpMap.getFreqAtSite(siteIndex);
//				int newStat = DATATYPE.getState(newChar);
//				double probOldToNew = 0;
//				double sum = 0;
//				for (int i = 0; i < charProb.length; i++) {
//					if(i != oldState){
//						sum += charProb[i];
//					}
//				}
				
				
				
//				newChar = getNextBase();
//			} while (siteType == SiteType.HIGH_VAR);	
		}
//		System.out.println(siteType +"\t"+ newChar +"\t"+ logq +"\t"+ 
//				Arrays.toString(srpMap.getFreqAtSite(siteIndex))	);
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
