package srp.operator.haplotypes;

import srp.evolution.OperationType;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

public class BasesReplaceOperator extends AbstractMultiOperator {


	public static final String OPERATOR_NAME = BasesReplaceOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.MULTI;
	

	
	public BasesReplaceOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
		super(haplotypeModel, length, mode);
		
	}


	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startAlignmentModelOperation();

		int hapIndex = getNextHapIndex();
		Haplotype haplotype = haplotypeModel.getHaplotype(hapIndex);
	    
		int sourceHapIndex = hapIndex;
	    do{
	    	sourceHapIndex = getNextHapIndex();
		} while(sourceHapIndex==hapIndex);

	    Haplotype sourceHap = haplotypeModel.getHaplotype(sourceHapIndex);
	    
	    int[] siteIndexs = generateUniqueSites(basesCount);

		for (int i : siteIndexs) {
			
			int oldState = haplotype.getState(i);
			int newState = sourceHap.getState(i);
			
			if(oldState != newState){
				char newChar = sourceHap.getChar(i);
				haplotype.setCharAt(i, newChar);
			}
		}
		
        // symmetrical move so return a zero hasting ratio
		haplotypeModel.setOperationRecord(OP, hapIndex, siteIndexs);
		
		haplotypeModel.endAlignmentModelOperation();
/*

AAAA
CCCC
GGGG
TTTT
ACGT	

to

AAAA
CCCC
GGGG
ACGT
ACGT

1/7
T -> A 1/4
T -> C 1/4 
T -> G 1/4
T -> T 1/4

1/256
A -> T 1/4
C -> T 1/4
G -> T 1/4
T -> T 1/4


*/
		
		
		
		
		return 0.0;//TODO This might be wrong!! should not equal for the reverse step
	}


	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

}
