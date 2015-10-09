package srp.operator.haplotypes;

import srp.evolution.OperationType;
import srp.evolution.haplotypes.HaplotypeModel;
import srp.evolution.shortreads.ShortReadMapping;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.math.MathUtils;

public abstract class AbstractHaplotypeOperator extends AbstractCoercableOperator {

	public static final DataType DATATYPE = Nucleotides.INSTANCE;
	public static final char[] DNA_CHARS = {'A','C','G','T'};
	public static final int DIMENSION = DNA_CHARS.length;
	
	public final int haplotypeCount;
	public final int haplotypeLength;
	protected HaplotypeModel haplotypeModel;
	protected ShortReadMapping srpMap;
	
	public AbstractHaplotypeOperator(HaplotypeModel haplotypeModel) {
		this(haplotypeModel, CoercionMode.COERCION_OFF);
//		this.haplotypeModel = haplotypeModel;
//		haplotypeCount = this.haplotypeModel.getHaplotypeCount();
//		haplotypeLength = this.haplotypeModel.getHaplotypeLength();
	}

	public AbstractHaplotypeOperator(HaplotypeModel haplotypeModel,
			CoercionMode mode) {
		super(mode);
		this.haplotypeModel = haplotypeModel;
		haplotypeCount = this.haplotypeModel.getHaplotypeCount();
		haplotypeLength = this.haplotypeModel.getHaplotypeLength();
		srpMap = this.haplotypeModel.getShortReadMapping();
	}

	public int getNextHapIndex() {
		return MathUtils.nextInt(haplotypeCount);
	}

	public int getNextSiteIndex() {
		return getNextSiteIndex(haplotypeLength);
	}
	public int getNextSiteIndex(int length){
		return MathUtils.nextInt(length);
	}
	
	public abstract OperationType getOperationType();

	public static char getNextDiffBase(int oldState) {
		int i = oldState;
		do {
			i = MathUtils.nextInt(DIMENSION);
		} while (i == oldState);
	
		return DNA_CHARS[i];
	}

	public static char getNextBase() {
		int i = MathUtils.nextInt(DIMENSION);
		return DNA_CHARS[i];
	}


	protected char transition(char oldChar){
		char newChar = oldChar;
		switch (oldChar){
		case 'A':
			newChar = 'G';
			break;
		case 'G':
			newChar = 'A';
			break;
		case 'C':
			newChar = 'T';
			break;
		case 'T':
			newChar = 'C';
			break;
		default:
			throw new IllegalArgumentException("Invalid char: "+oldChar);
		}
		
		
		return newChar;
	}
	
	char[] tr_AG = {'C', 'T'};
	char[] tr_CT = {'A', 'G'};
	
	protected char transversion(char oldChar) {
		int i = MathUtils.nextInt(2);
		char newChar = oldChar;
		switch (oldChar){
		case 'A':
		case 'G':
			newChar = tr_AG[i];
			break;
//		case 'G':
//			newChar = 'A';
//			break;
		case 'C':
		case 'T':
			newChar = tr_CT[i];
			break;
//		case 'T':
//			newChar = 'C';
//			break;
		default:
			throw new IllegalArgumentException("Invalid char: "+oldChar);
		}
		
		
		return newChar;
	}

	protected boolean checkUnique(int i) {
		
		char baseChar = haplotypeModel.getHaplotype(0).getChar(i);
		for (int j = 1; j < haplotypeCount; j++) {
			char charJ = haplotypeModel.getHaplotype(j).getChar(i);
			if(charJ != baseChar){
				return true;
			}
		}
		return false;
	}


	
	
}