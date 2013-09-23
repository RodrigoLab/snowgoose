package srp.haplotypes.operator;

import java.util.Arrays;

import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.SwapInfo;
import dr.evolution.datatype.Nucleotides;
import dr.inference.model.Parameter;
import dr.inference.operators.OperatorFailedException;

public class SingleBaseFrequencyOperator extends AbstractSingleBaseOperator {

	public final static String OPERATOR_NAME = SingleBaseFrequencyOperator.class.getSimpleName();
	private static final int NEW_CHAR_INDEX = SwapInfo.SWAP_BASE_NEW_CHAR_INDEX;
	private static final int OLD_CHAR_INDEX = SwapInfo.SWAP_BASE_OLD_CHAR_INDEX;
	
	private Parameter frequency;

	public SingleBaseFrequencyOperator(HaplotypeModel haplotypeModel,
			Parameter freqs) {

		super(haplotypeModel);
		this.frequency = freqs;
	}

	/*
    *
    * @return the log-transformed hastings ratio
    */
	@Override
	public double doOperation() throws OperatorFailedException {

		haplotypeModel.startHaplotypeOperation();
		
		int[] posChar = haplotypeModel.getNextBaseFrequency(frequency);

		int[] swapRecord = haplotypeModel.swapHaplotypeSingleBase(OP, posChar);
		
		haplotypeModel.endHaplotypeOperation();
		
		if(swapRecord[OLD_CHAR_INDEX]==swapRecord[NEW_CHAR_INDEX]){
			return 0.0;
		}
		else{
			int newChar = Nucleotides.NUCLEOTIDE_STATES[swapRecord[NEW_CHAR_INDEX]];
			int oldChar = Nucleotides.NUCLEOTIDE_STATES[swapRecord[OLD_CHAR_INDEX]];
			double oldProb = frequency.getParameterValue(oldChar);
			double newProb = frequency.getParameterValue(newChar);
			double logq = Math.log(oldProb/newProb);
//			System.out.println(newChar +"\t"+ oldChar +"\t"+ logq+"\t"+ (oldProb/newProb) +"\t"+ Arrays.toString(swapRecord) +"\t"+  Arrays.toString(frequency.getParameterValues()));
			

			return logq;
		}
	}

	public static final int NUCLEOTIDE_STATES[] = {
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 0-15
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 16-31
	//                                          -
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 32-47
	//                                                ?
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,16,	// 48-63
	//	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
		17, 0,11, 1,12,16,16, 2,13,16,16,10,16, 7,15,16,	// 64-79
	//	 p  q  R  S  T  U  V  W  x  Y  z
		16,16, 5, 9, 3, 3,14, 8,16, 6,16,17,17,17,17,17,	// 80-95
	//	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
		17, 0,11, 1,12,16,16, 2,13,16,16,10,16, 7,15,16,	// 96-111
	//	 p  q  R  S  T  U  V  W  x  Y  z
		16,16, 5, 9, 3, 3,14, 8,16, 6,16,17,17,17,17,17		// 112-127
	};
	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

}
