package srp.operator.haplotypes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import srp.evolution.OperationType;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class ParsimonyInformativeBasesMultiOperator extends AbstractMultiOperator {


	public static final String OPERATOR_NAME = ParsimonyInformativeBasesMultiOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.MULTI;
	

	
	public ParsimonyInformativeBasesMultiOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
		super(haplotypeModel, length, mode);
		
	}


	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startAlignmentModelOperation();

		int hapIndex = getNextHapIndex();
		Haplotype haplotype = haplotypeModel.getHaplotype(hapIndex);
//	    
//	    
		double xNewToOld = 1.0;// / Math.pow(4, basesCount);
		double xOldToNew = 1.0;// / (haplotypeCount-1);
		
		int[] siteIndexs = new int[basesCount];
//		int basesCountMinusOne = basesCount-1;
		int count = 0;
		siteIndexs[count] = MathUtils.nextInt(haplotypeLength);  
		int site = siteIndexs[count];
		
		count++;
		
		do{
			
			site++;
			if(site == haplotypeLength){
				site = 0;
			}
				
			int[] parsimonyInfo = checkParsimony(site);
			if(parsimonyInfo!=null){

				siteIndexs[count] = site;
				count++;
				char oldChar = haplotype.getChar(site);
				char newChar;// = oldChar;
				int oldState = haplotype.getState(site);
				
				int s = oldState;
				do {
					s = MathUtils.nextInt(DIMENSION);
					newChar = DNA_CHARS[s];
				} while ((s == oldState) || (parsimonyInfo[newChar]==0) );
				xNewToOld *= 1.0 / (haplotypeCount-parsimonyInfo[oldChar]);
				xOldToNew *= 1.0 / (haplotypeCount-parsimonyInfo[newChar]);
				
				
//					System.out.println( (oldChar != newChar) +"\t"+ (parsimonyInfo[oldChar]>0) +"\t"+ (parsimonyInfo[newChar]>0) );
//					System.out.println( oldChar +"\t"+ newChar +"\t"+ (parsimonyInfo[oldChar]) +"\t"+ (parsimonyInfo[newChar]) );
//					System.out.println(xNewToOld +"\t"+  xOldToNew +"\t"+ (xNewToOld / xOldToNew) );
//					System.out.println(xNewToOld +"\t"+  xOldToNew +"\t"+ (xNewToOld / xOldToNew) +"\t"+ haplotypeCount +"\t"+ Arrays.toString(haplotypeModel.getSitePattern(i))+"\n");
				
				haplotype.setCharAt(site, newChar);
			}
//				else{
//					System.out.println("=====\t"+Arrays.toString(haplotypeModel
//							.getSitePattern(i)));
//				}
			if (count == basesCount){
				break;
			}
	//	        System.out.println(hapIndex +"\t"+ i +"\t from "+oldState +" to new "+ newChar);
			
		}
		while(count<basesCount);
		
		HashSet<Integer> hs = new HashSet<Integer>(  );
		for (int i = 0; i < siteIndexs.length; i++) {
			hs.add(siteIndexs[i]);
		}
		if(hs.size() < (basesCount-1) ){
			System.out.println("REPEAT:"+Arrays.toString(siteIndexs));
		}
		
//		System.out.println("SITE:"+Arrays.toString(siteIndexs) +"\t"+ count);
        // symmetrical move so return a zero hasting ratio
		haplotypeModel.setOperationRecord(OP, hapIndex, siteIndexs);
	
		haplotypeModel.endAlignmentModelOperation();

//		double xNewToOld = 1.0 / Math.pow(4, basesCount);
//		double xOldToNew = 1.0 / (haplotypeCount-1);
		
		double logq = Math.log( xNewToOld ) - Math.log( xOldToNew  );
//		System.out.println("logQ:\t"+xNewToOld +"\t"+  xOldToNew +"\t"+ (xNewToOld / xOldToNew) +"\t"+ logq);
		return logq;//TODO This might be wrong!! should not equal for the reverse step
		
	}


	private int[] checkParsimony(int i) {
		int[] charCount = new int['U'];
//		char baseChar = haplotypeModel.getHaplotype(0).getChar(i);
		for (int h = 0; h < haplotypeCount; h++) {
			char charH = haplotypeModel.getHaplotype(h).getChar(i);
			charCount[charH]++;
		}
		int parInfoCount = 0;
		if (charCount['A']>1){
			parInfoCount++;
		}
		if (charCount['C']>1){
			parInfoCount++;
		}
		if (charCount['G']>1){
			parInfoCount++;
		}
		if (charCount['T']>1){
			parInfoCount++;
		}
		if(parInfoCount>1){
//			System.out.println(Arrays.toString(charCount));
			return charCount;
		}
		return null;
	}


	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

}
