package srp.operator.haplotypes;

import java.util.Arrays;

import srp.evolution.OperationType;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class PopulateSingletonOperator extends AbstractMultiOperator {


	public static final String OPERATOR_NAME = PopulateSingletonOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.MULTI;
	

	
	public PopulateSingletonOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
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
		do{
//			System.out.println("loop:"+count);
			int[] tempSiteIndexs = generateUniqueSites(basesCount*10);
			for (int i : tempSiteIndexs) {
				
				int[] parsimonyInfo = checkSingleton(i);
				if(parsimonyInfo!=null){

					siteIndexs[count] = i;
					count++;
					char oldChar = haplotype.getChar(i);
					char newChar;// = oldChar;
					int oldState = haplotype.getState(i);
					
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
					
					haplotype.setCharAt(i, newChar);
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
		}
		while(count<basesCount);
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


	private int[] checkSingleton(int i) {
		int[] charCount = new int['U'];
//		char baseChar = haplotypeModel.getHaplotype(0).getChar(i);
		for (int h = 0; h < haplotypeCount; h++) {
			char charH = haplotypeModel.getHaplotype(h).getChar(i);
			charCount[charH]++;
		}
		int parInfoCount = 0;
		if (charCount['A']==1){
			parInfoCount++;
		}
		if (charCount['C']==1){
			parInfoCount++;
		}
		if (charCount['G']==1){
			parInfoCount++;
		}
		if (charCount['T']==1){
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
