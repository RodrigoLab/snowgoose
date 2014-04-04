package srp.likelihood.spectrum;

import java.util.Set;

import srp.evolution.OperationType;
import srp.evolution.shortreads.ShortReadMapping;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import srp.likelihood.AbstractShortReadsLikelihood;
import srp.likelihood.stateLikelihood.BetaMeanStateLikelihood;
import srp.likelihood.stateLikelihood.BetaModeStateLikelihood;
import srp.likelihood.stateLikelihood.ChisqStateLikelihood;
import srp.likelihood.stateLikelihood.GTestStateLikelihood;
import srp.likelihood.stateLikelihood.ProbabilityStateLikelihood;
import srp.likelihood.stateLikelihood.StateLikelihood;

import com.carrotsearch.hppc.BitSet;

public abstract class AbstractShortReadsSpectrumLikelihood extends AbstractShortReadsLikelihood {

	private static final long serialVersionUID = 2079474866153379297L;

	protected StateLikelihood stateLikelihood;
	protected SpectrumAlignmentModel spectrumModel;
	


	protected double[] eachSrpLogLikelihood;
	protected double[] storedEachSrpLogLikelihood;

	protected double[] sumScaledSrpLogLikelihood;
	protected double[] storedSumSrpLogLikelihood;
	
	
	protected double[] spectrumLogLikelihood;
	protected double[] storedSpectrumLogLikelihood;
	
	protected double[] scaledSpectrumLogLikelihood;
	protected double[] storedScaledSpectrumLogLikelihood;
	
	
	public AbstractShortReadsSpectrumLikelihood(String name) {
		super(name);
	}


	protected void setDistType(DistType distType) {
		// this.distType = distType;
		try {

			switch (distType) {
			case flat:
				stateLikelihood = new ProbabilityStateLikelihood();
				break;
			case betaMean:
				stateLikelihood = new BetaMeanStateLikelihood();
				break;
			case betaMode:
				stateLikelihood = new BetaModeStateLikelihood();
				break;
			case gTest:
				stateLikelihood = new GTestStateLikelihood();
				break;
			case chisq:
				stateLikelihood = new ChisqStateLikelihood();
				break;
			default:
				throw new IllegalArgumentException("Invalid distType: "
						+ distType);
			}
			likelihoodKnown = false;

		} catch (IllegalArgumentException e) {
			System.err.println("Invalid distribution type " + distType);
			e.printStackTrace();
			System.exit(-1);
		}
	}

	protected void storeIJ(int i, int j) {
		
		int offset = i*spectrumCount+j;
		storedSpectrumLogLikelihood[offset] = spectrumLogLikelihood[offset];
		storedScaledSpectrumLogLikelihood[offset] = scaledSpectrumLogLikelihood[offset];
	}

	protected void storeI(int i) {
		storedEachSrpLogLikelihood[i] = eachSrpLogLikelihood[i];
		storedSumSrpLogLikelihood[i] = sumScaledSrpLogLikelihood[i];
	}
	
	protected void restoreIJ(int i, int j) {
		int offset = i*spectrumCount+j;
		spectrumLogLikelihood[offset] = storedSpectrumLogLikelihood[offset];
		scaledSpectrumLogLikelihood[offset] = storedScaledSpectrumLogLikelihood[offset];

//		spectrumLogLikelihood2D[i][j] = storedSpectrumLogLikelihood2D[i][j];
//		scaledSpectrumLogLikelihood2D[i][j] = storedScaledSpectrumLogLikelihood2D[i][j];
	}

	protected void restoreI(int i) {
		eachSrpLogLikelihood[i] = storedEachSrpLogLikelihood[i];
		sumScaledSrpLogLikelihood[i] = storedSumSrpLogLikelihood[i]; 
		
	}



	protected void storeEverything() {
	
		System.arraycopy(eachSrpLogLikelihood, 0, storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		System.arraycopy(sumScaledSrpLogLikelihood, 0, storedSumSrpLogLikelihood, 0, sumScaledSrpLogLikelihood.length);
		System.arraycopy(spectrumLogLikelihood,0, storedSpectrumLogLikelihood, 0, spectrumLogLikelihood.length);
		System.arraycopy(scaledSpectrumLogLikelihood,0, storedScaledSpectrumLogLikelihood, 0, scaledSpectrumLogLikelihood.length);
	}

	protected void restoreEverything(){
		
		System.arraycopy(storedEachSrpLogLikelihood, 0, eachSrpLogLikelihood, 0, eachSrpLogLikelihood.length);
		System.arraycopy(storedSumSrpLogLikelihood, 0, sumScaledSrpLogLikelihood, 0, sumScaledSrpLogLikelihood.length);
		System.arraycopy(storedSpectrumLogLikelihood, 0, spectrumLogLikelihood, 0, storedSpectrumLogLikelihood.length);
		System.arraycopy(storedScaledSpectrumLogLikelihood, 0, scaledSpectrumLogLikelihood, 0, storedScaledSpectrumLogLikelihood.length);
	}
	

	@Override
	public void makeDirty() {
		// System.err.println("make dirty");
		spectrumModel.resetOperation();
		likelihoodKnown = false;
	
	}
	
	public OperationType getOperation() {
		return spectrumModel.getOperation();
	}

	
	public enum DistType {
		betaMean(0), betaMode(1), gTest(2), chisq(3), flat(9), ;
		int code;

		private DistType(int code) {
			this.code = code;
		}
		// DistType.valueOf(codeString);
	}

}