package srp.likelihood.spectrum;

import java.util.ArrayList;
import java.util.Set;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.carrotsearch.hppc.BitSet;

import srp.likelihood.stateLikelihood.BetaMeanStateLikelihood;
import srp.likelihood.stateLikelihood.BetaModeStateLikelihood;
import srp.likelihood.stateLikelihood.ChisqStateLikelihood;
import srp.likelihood.stateLikelihood.GTestStateLikelihood;
import srp.likelihood.stateLikelihood.ProbabilityStateLikelihood;
import srp.likelihood.stateLikelihood.StateLikelihood;
import srp.shortreads.ShortReadMapping;
import srp.spectrum.SpectrumAlignmentModel;
import srp.spectrum.SpectrumOperation;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;

public abstract class AbstractShortReadsSpectrumLikelihood extends
		AbstractModelLikelihood {

	private static final long serialVersionUID = 2079474866153379297L;


	protected double logLikelihood;
	protected double storedLogLikelihood;
	
	
	protected boolean likelihoodKnown;
	
	protected int spectrumLength;
	protected int spectrumCount;
	protected StateLikelihood stateLikelihood;
	protected MultiType multiType;

	protected SpectrumAlignmentModel spectrumModel;
	protected ShortReadMapping srpMap;
	
	protected int[][] mapToSrpArray;

	
	protected boolean[] srpSwitch;
	protected Set<Integer> allSrpPos;
	protected BitSet bitSet;
	
	public AbstractShortReadsSpectrumLikelihood(String name) {
		super(name);
	}

	
	protected double getStoredLogLikelihood() {
		return storedLogLikelihood;
	}
	
	protected abstract double calculateLogLikelihood();


	protected void recalculateArray(int[] siteIndexs) {
		for (int s : siteIndexs) {
			for (int i : mapToSrpArray[s]){
				srpSwitch[i] = true;
			}
		}
	}

	protected void recalculateHashSet(int[] siteIndexs) {
		allSrpPos.clear();
		for (int s : siteIndexs) {
			ArrayList<Integer> mapToSrp = srpMap.getMapToSrp(s);
			allSrpPos.addAll(mapToSrp);
		}
	}

	protected void recalculateBitSet(int[] siteIndexs) {
		bitSet.clear();
//		BitSet bitSet = new BitSet(srpCount);
		for (int s : siteIndexs) {
			BitSet tempSet = srpMap.getBitSet(s);
			bitSet.or(tempSet);
		}
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

	public SpectrumOperation getOperation() {
		return spectrumModel.getSpectrumOperation();
	}


	public void setMultiType(MultiType type) {
		this.multiType = type;
	}


	@Override
	public double getLogLikelihood() {
	
		if (!likelihoodKnown) {
			logLikelihood = calculateLogLikelihood();
			likelihoodKnown = true;
		}
	
		return logLikelihood;
	
	}


	@Override
	public Model getModel() {
		return this;
	
	}


	@Override
	public void makeDirty() {
		// System.err.println("make dirty");
		spectrumModel.resetSpectrumOperation();
		likelihoodKnown = false;
	
	}

	@Override
	public Element createElement(Document d) {
		throw new RuntimeException("Not implemented yet!");
	}

	protected enum MultiType {
		Array, Hash, All, BitSet, ;

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