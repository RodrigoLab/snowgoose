package srp.operator.haplotypes;

import srp.evolution.OperationType;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

public class BasesConsecutiveOperator extends AbstractMultiOperator {


	public static final String OPERATOR_NAME = BasesConsecutiveOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.MULTI;
	

	private int maxLength;
	
	public BasesConsecutiveOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
		super(haplotypeModel, length, mode);
		maxLength = haplotypeLength-basesCount;
	}


	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startAlignmentModelOperation();

		int hapIndex = getNextHapIndex();
		Haplotype haplotype = haplotypeModel.getHaplotype(hapIndex);
	    maxLength = haplotypeLength - basesCount;
	    int siteStart = getNextSiteIndex(maxLength);
		int[] siteIndexs = new int[basesCount];
		for (int i = 0; i < siteIndexs.length; i++) {
			siteIndexs[i] = siteStart + i; 
		}
		double logq = 0;
		for (int i : siteIndexs) {
			
//			SpectraParameter spectra = spectrum.getSpectra(siteIndexs[i]);
//			swapFrequency(spectra);
//			int oldState = haplotype.getState(i);
//			char newChar = getNextDiffBase(oldState);
//			newChar = getNextBase();//TODO: which one?
			
			char newChar = srpMap.nextBaseFreqAt(i);
			char oldChar = haplotype.getChar(i);

			logq += srpMap.calculateLogqOldNewChar(i, oldChar, newChar);
			
			haplotype.setCharAt(i, newChar);

//	        System.out.println(hapIndex +"\t"+ i +"\t from "+oldState +" to new "+ newChar);
		}
        // symmetrical move so return a zero hasting ratio
		haplotypeModel.setOperationRecord(OP, hapIndex, siteIndexs);
	
		haplotypeModel.endAlignmentModelOperation();

		return logq;
	}


	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

}
