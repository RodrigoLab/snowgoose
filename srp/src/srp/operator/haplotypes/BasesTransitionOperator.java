package srp.operator.haplotypes;

import srp.evolution.OperationType;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

public class BasesTransitionOperator extends AbstractMultiOperator {


	public static final String OPERATOR_NAME = BasesTransitionOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.MULTI;
	private int maxLength;
	

	
	public BasesTransitionOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
		super(haplotypeModel, length, mode);
		maxLength = haplotypeLength-basesCount;
	}


	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startAlignmentModelOperation();

		int hapIndex = getNextHapIndex();
		Haplotype haplotype = haplotypeModel.getHaplotype(hapIndex);
//	    
	    int siteStart = getNextSiteIndex(maxLength);
		int[] siteIndexs = new int[basesCount];
		for (int i = 0; i < siteIndexs.length; i++) {
			siteIndexs[i] = siteStart + i; 
		}

		for (int i : siteIndexs) {
			
//			SpectraParameter spectra = spectrum.getSpectra(siteIndexs[i]);
//			swapFrequency(spectra);
			
			char oldState = haplotype.getChar(i);
			char newChar = transition(oldState);
			
			haplotype.setCharAt(i, newChar);
			
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
