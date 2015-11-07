package srp.operator.haplotypes;

import srp.evolution.OperationType;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

public class BasesDataMultiOperator extends AbstractMultiOperator {


	public static final String OPERATOR_NAME = BasesDataMultiOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.MULTI;
	

	
	public BasesDataMultiOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
		super(haplotypeModel, length, mode);
		
	}


	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startAlignmentModelOperation();

		int hapIndex = getNextHapIndex();
		Haplotype haplotype = haplotypeModel.getHaplotype(hapIndex);
//	    
//	    
		int[] siteIndexs = generateUniqueSites(basesCount);

		for (int i : siteIndexs) {
			
//			boolean isFix;
//			do {
//				i= getNextSiteIndex();
//				isFix = srpMap.isFixSite(i);
////				System.out.println("DO: "+siteIndex +"\t"+ isFix);
//			} while (isFix);
			
//			SpectraParameter spectra = spectrum.getSpectra(siteIndexs[i]);
//			swapFrequency(spectra);
//			int oldState = haplotype.getState(i);
//			char newChar = getNextDiffBase(oldState);
			
//			newChar = getNextBase();
			char newChar = srpMap.nextBaseFreqAt(i);
			haplotype.setCharAt(i, newChar);
//	        System.out.println(hapIndex +"\t"+ i +"\t from "+oldState +" to new "+ newChar);
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
