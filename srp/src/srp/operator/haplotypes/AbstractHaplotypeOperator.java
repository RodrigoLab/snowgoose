package srp.operator.haplotypes;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;

import srp.evolution.OperationType;
import srp.haplotypes.HaplotypeModel;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

public abstract class AbstractHaplotypeOperator extends AbstractCoercableOperator {

	public static final DataType DATATYPE = Nucleotides.INSTANCE;
	public static final char[] DNA_CHARS = {'A','C','G','T'};
	public static final int DIMENSION = DNA_CHARS.length;
	
	public final int haplotypeCount;
	public final int haplotypeLength;
	protected HaplotypeModel haplotypeModel;

	public AbstractHaplotypeOperator(HaplotypeModel haplotypeModel) {
		super(CoercionMode.COERCION_OFF);
		this.haplotypeModel = haplotypeModel;
		haplotypeCount = this.haplotypeModel.getHaplotypeCount();
		haplotypeLength = this.haplotypeModel.getHaplotypeLength();
	}

	public AbstractHaplotypeOperator(HaplotypeModel haplotypeModel,
			CoercionMode mode) {
		super(mode);
		this.haplotypeModel = haplotypeModel;
		haplotypeCount = this.haplotypeModel.getHaplotypeCount();
		haplotypeLength = this.haplotypeModel.getHaplotypeLength();
	}

	public int getNextHapIndex() {
		return MathUtils.nextInt(haplotypeCount);
	}

	public int getNextSiteIndex() {
		return MathUtils.nextInt(haplotypeLength);
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

}