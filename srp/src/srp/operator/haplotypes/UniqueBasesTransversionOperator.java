package srp.operator.haplotypes;

import java.util.Arrays;

import srp.evolution.OperationType;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;
import srp.operator.haplotypes.BasesTransitionOperator;

public class UniqueBasesTransversionOperator extends AbstractMultiOperator {


	public static final String OPERATOR_NAME = UniqueBasesTransversionOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.MULTI;
	

	
	public UniqueBasesTransversionOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
		super(haplotypeModel, length, mode);
		
	}


@Override
		public double doOperation() throws OperatorFailedException {
	
			haplotypeModel.startAlignmentModelOperation();
	
			int hapIndex = getNextHapIndex();
			Haplotype haplotype = haplotypeModel.getHaplotype(hapIndex);
			
			int[] siteIndexs = new int[basesCount];
//			int basesCountMinusOne = basesCount-1;
			int count = 0;
			siteIndexs[count] = MathUtils.nextInt(haplotypeLength);  
			int site = siteIndexs[count];
			
//			count++;
			
			do{
				
				site++;
				if(site == haplotypeLength){
					site = 0;
				}
					
					
				boolean isNotUnique = !checkUnique(site);
				if(isNotUnique){

//					System.out.println(Arrays.toString(haplotypeModel.getSitePattern(i)));
					siteIndexs[count] = site;
					count++;
//					int oldState = haplotype.getState(site);
//					char newChar = getNextDiffBase(oldState);
////					newChar = getNextBase();
//					haplotype.setCharAt(site, newChar);
					char oldChar = haplotype.getChar(site);
					char newChar = transversion(oldChar);
					
					haplotype.setCharAt(site, newChar);
					
				}
//				else{
//					System.out.println("=====\t"+Arrays.toString(haplotypeModel.getSitePattern(i)));
//				}
				if (count == basesCount){
					break;
				}
	//	        System.out.println(hapIndex +"\t"+ i +"\t from "+oldState +" to new "+ newChar);
				
			}
			while(count<basesCount);
			
	        // symmetrical move so return a zero hasting ratio
			haplotypeModel.setOperationRecord(OP, hapIndex, siteIndexs);
		
			haplotypeModel.endAlignmentModelOperation();
	
			return 0.0;
		}


	//	@Override
	public double doOperation2() throws OperatorFailedException {

		haplotypeModel.startAlignmentModelOperation();

		int hapIndex = getNextHapIndex();
		Haplotype haplotype = haplotypeModel.getHaplotype(hapIndex);
		
		int[] siteIndexs = new int[basesCount];
		int count = 0;
		do{
//			System.out.println("loop:"+count);
			int[] tempSiteIndexs = generateUniqueSites(basesCount*10);
			for (int i : tempSiteIndexs) {
				
				boolean isNotUnique = checkUnique(i);
				if(isNotUnique){

//					System.out.println(Arrays.toString(haplotypeModel
//							.getSitePattern(i)));
					siteIndexs[count] = i;
					count++;
					int oldState = haplotype.getState(i);
					char newChar = getNextDiffBase(oldState);
//					newChar = getNextBase();
					haplotype.setCharAt(i, newChar);
					break;
				}
			}
		}
		while(count<1);
		
		
		for (int i = siteIndexs[0]; i < haplotypeLength; i++) {
			
		
			boolean isNotUnique = checkUnique(i);
			if(isNotUnique){

//				System.out.println(Arrays.toString(haplotypeModel
//						.getSitePattern(i)));
				siteIndexs[count] = i;
				count++;
				int oldState = haplotype.getState(i);
				char newChar = getNextDiffBase(oldState);
//				newChar = getNextBase();
				haplotype.setCharAt(i, newChar);
			}
			if (count == basesCount){
				break;
			}
		}
		
        // symmetrical move so return a zero hasting ratio
		haplotypeModel.setOperationRecord(OP, hapIndex, siteIndexs);
	
		haplotypeModel.endAlignmentModelOperation();

		return 0.0;
	}


	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

}
