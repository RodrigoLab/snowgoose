package srp.operator.haplotypes;

import srp.evolution.OperationType;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;

import java.util.Arrays;

import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public class BasesDataMultiOperator extends AbstractMultiOperator {


	public  String OPERATOR_NAME;// = BasesDataMultiOperator.class.getSimpleName();
	public static final OperationType OP = OperationType.MULTI;
	private int type = 3;
	

	
	public BasesDataMultiOperator(HaplotypeModel haplotypeModel, int length, CoercionMode mode) {
		super(haplotypeModel, length, mode);
	}
	
	public BasesDataMultiOperator(HaplotypeModel haplotypeModel, double pct, CoercionMode mode) {
		super(haplotypeModel, (int) (haplotypeModel.getHaplotypeLength()*
				((pct<1)? pct : 0.01) ), mode);
	}

	public BasesDataMultiOperator(HaplotypeModel haplotypeModel, double pct, CoercionMode mode, 
			int type) {
		super(haplotypeModel, (int) (haplotypeModel.getHaplotypeLength()*
				((pct<1)? pct : 0.01) ), mode);
		this.type = type;
		OPERATOR_NAME = BasesDataMultiOperator.class.getSimpleName()+"_"+this.type;
	}
	
	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startAlignmentModelOperation();

		int hapIndex = getNextHapIndex();
		Haplotype haplotype = haplotypeModel.getHaplotype(hapIndex);
//	    
//	    
//		int[] siteIndexs = generateUniqueSites(basesCount);
		int[] siteIndexs;
		switch (type) {
		case 1:
			siteIndexs = randomSampleSites(basesCount, srpMap.highVarSiteArray);
			break;
		case 2:
			siteIndexs = randomSampleSites(basesCount, srpMap.lowVarSiteArray);
			break;
		case 3:
			 siteIndexs = generateUniqueSites(basesCount);
			break;
		default:
			siteIndexs = generateUniqueSites(basesCount);
			break;
		}
		
		
		
//		System.out.println(Arrays.toString(srpMap.highVarSiteArray));
//		System.out.println("\t"+Arrays.toString(siteIndexs));

		double logq = 0;
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
