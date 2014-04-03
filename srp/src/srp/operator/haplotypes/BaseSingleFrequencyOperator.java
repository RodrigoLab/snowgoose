package srp.operator.haplotypes;

import srp.haplotypes.HaplotypeModel;
import srp.haplotypes.SwapInfo;
import dr.evolution.datatype.Nucleotides;
import dr.inference.model.Parameter;
import dr.inference.operators.OperatorFailedException;

public class BaseSingleFrequencyOperator extends AbstractBaseSingleOperator {

	
	public final static String OPERATOR_NAME = BaseSingleFrequencyOperator.class.getSimpleName();
	
	private static final int NUCLEOTIDE_STATES[] = Nucleotides.NUCLEOTIDE_STATES;
	private static final int NEW_CHAR_INDEX = SwapInfo.SWAPBASE_NEW_CHAR_INDEX;
	private static final int OLD_CHAR_INDEX = SwapInfo.SWAPBASE_OLD_CHAR_INDEX;
	
	private Parameter frequency;

	public BaseSingleFrequencyOperator(HaplotypeModel haplotypeModel,
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
		
		double logq = 0.0;
		if(swapRecord[OLD_CHAR_INDEX]!=swapRecord[NEW_CHAR_INDEX]){

//			int newChar = NUCLEOTIDE_STATES[swapRecord[NEW_CHAR_INDEX]];
//			int oldChar = NUCLEOTIDE_STATES[swapRecord[OLD_CHAR_INDEX]];
			logq = haplotypeModel.getLogqFrequency(swapRecord[OLD_CHAR_INDEX], swapRecord[NEW_CHAR_INDEX]);

		}
		return logq;
	}

	
	@Override
	public String getOperatorName() {
	
		return OPERATOR_NAME;
	}

}
