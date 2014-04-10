package srp.haplotypes;


import java.util.List;

import jebl.evolution.sequences.Sequence;
import srp.dr.ext.SeqGenExt;
import srp.dr.ext.TreeLikelihoodExt;
import srp.evolution.OperationType;
import srp.evolution.shortreads.AlignmentMapping;
import srp.evolution.spectrum.Spectrum;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.util.NumberFormatter;

public class HaplotypeModel extends AbstractHaplotypeModel  {
	
	private static final long serialVersionUID = -5057514703825711955L;
	
	private static final String MODEL_NAME = "HaplotypeModel";
	private static final int NUCLEOTIDE_STATES[] = Nucleotides.NUCLEOTIDE_STATES;
	private static final int STATE_COUTN = DATA_TYPE.getStateCount();

	private boolean DEBUG = false;


	private void initHaplotypes() {
		for (int i = 0; i < sequenceCount; i++) {
			Taxon taxon = new Taxon(TAXON_PREFIX + i);
			Haplotype haplotype = new Haplotype(taxon, sequenceLength);
			setHaplotype(i, haplotype);
		}

	}
	
	public HaplotypeModel(String[] sequences){
		this(AlignmentUtils.createAlignment(sequences));
	}

	
	public HaplotypeModel(int hapCount, int hapLength) {
		super(MODEL_NAME, hapCount, hapLength);
		initHaplotypes();
	}
	
	public HaplotypeModel(Alignment trueAlignment) {
		this(trueAlignment.getSequenceCount(), trueAlignment.getSiteCount());

		for (int i = 0; i < trueAlignment.getSequenceCount(); i++) {
			Haplotype haplotype = new Haplotype(trueAlignment.getSequence(i));
			setHaplotype(i, haplotype);
		}
	}
	
	
	public int calculateSPS(){//TODO test this
		return SPSDist.calculeteSPS(this, this);
	}

	@Override
	public String toString(){

        NumberFormatter formatter = new NumberFormatter(6);

        StringBuilder buffer = new StringBuilder();

//	        boolean countStatistics = !(dataType instanceof Codons) && !(dataType instanceof GeneralDataType);

//        if (countStatistics) {
//            buffer.append("Site count = ").append(getSiteCount()).append("\n");
//            buffer.append("Invariant sites = ").append(getInvariantCount()).append("\n");
//            buffer.append("Singleton sites = ").append(getSingletonCount()).append("\n");
//            buffer.append("Parsimony informative sites = ").append(getInformativeCount()).append("\n");
//            buffer.append("Unique site patterns = ").append(getUniquePatternCount()).append("\n\n");
//        }
        for (int i = 0; i < getHaplotypeCount(); i++) {
            String name = formatter.formatToFieldWidth(getTaxon(i).getId(), 10);
            buffer.append(">" + name + "\n");
            buffer.append(getAlignedSequenceString(i) + "\n");
        }

        return buffer.toString();
    
		
	}

	@Override
	protected void storeState() {

		OperationType operation = operationRecord.getOperation();
		int spectrumIndex;
		int siteIndex;
		int[] siteIndexs;
		Haplotype haplotype;
		switch (operation) {
		case NONE:
			if (DEBUG) {
				System.out.println("StoreState in SpectrumAlignment:\t"
						+ operation);
			}
			break;
		case FULL:
			if (DEBUG) {
				System.out.println("StoreState in SpectrumAlignment:\t"
						+ operation);
			}
			for (int i = 0; i < getHaplotypeCount(); i++) {
				haplotype = getHaplotype(i);
				haplotype.storeState();
			}
			break;

		case SINGLE:
			spectrumIndex = operationRecord.getSpectrumIndex();
			siteIndex = operationRecord.getSingleIndex();
			haplotype = getHaplotype(spectrumIndex);
			haplotype.storeState(siteIndex);
			break;

		case MULTI:
			spectrumIndex = operationRecord.getSpectrumIndex();
			siteIndexs = operationRecord.getAllSiteIndexs();
			haplotype = getHaplotype(spectrumIndex);
			for (int s = 0; s < siteIndexs.length; s++) {
				haplotype.storeState(s);
			}
			break;

		case COLUMN:
			siteIndex = operationRecord.getSingleIndex();
			for (int i = 0; i < getHaplotypeCount(); i++) {
				haplotype = getHaplotype(i);
				haplotype.storeState(siteIndex);
			}
			break;

		case RECOMBINATION:
			int[] twoSpectrums = operationRecord.getRecombinationSpectrumIndex();
			int[] twoPositions = operationRecord.getRecombinationPositionIndex();
			for (int i : twoSpectrums) {
				haplotype = getHaplotype(i);
				for (int s = twoPositions[0]; s < twoPositions[1]; s++) {
					haplotype.storeState(s);
				}
			}
			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: "
					+ operation + "\tin"
					+ SpectrumAlignmentModel.class.getSimpleName());

		}
	}

	@Override
	protected void restoreState() {
	
		OperationType operation = operationRecord.getOperation();
		int hapIndex;
		int siteIndex;
		Haplotype haplotype;
		int[] siteIndexs;

//		System.err.println("zzzzRestore SpectrumAlignment: "+operation);
		switch (operation) {
		
		case NONE:
			if(DEBUG){
				System.out.println("RestoreState in SpectrumAlignment:\t"+operation);
			}
			break;
		case FULL:
			if(DEBUG){
				System.out.println("RestoreState in SpectrumAlignment:\t"+operation);
			}
			for (int i = 0; i < getHaplotypeCount(); i++) {
				haplotype = getHaplotype(i);
				haplotype.restoreState();
			}
			break;

		case SINGLE:
			hapIndex = operationRecord.getSpectrumIndex();
			siteIndex = operationRecord.getSingleIndex();
			haplotype = getHaplotype(hapIndex);
			haplotype.restoreState(siteIndex);
			break;
		
		case MULTI:
			hapIndex = operationRecord.getSpectrumIndex();
			siteIndexs = operationRecord.getAllSiteIndexs();
			haplotype = getHaplotype(hapIndex);
			for (int s: siteIndexs){
				haplotype.restoreState(s);
			}
			break;

		case COLUMN:
	//			spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			siteIndex = operationRecord.getSingleIndex();
			for (int i = 0; i < getHaplotypeCount(); i++) {
				haplotype = getHaplotype(i);
				haplotype.restoreState(siteIndex);
			}
			break;
			
		case RECOMBINATION:
//			System.err.println("Restore alignment recombination");
			int[] twoSpectrums = operationRecord.getRecombinationSpectrumIndex();
			int[] twoPositions = operationRecord.getRecombinationPositionIndex();
			for (int i : twoSpectrums) {
				haplotype = getHaplotype(i);
				for (int s = twoPositions[0]; s < twoPositions[1]; s++) {
					haplotype.restoreState(s);
				}
			}
			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: " + operation);

		}
	}
	

	public void simulateSequence(TreeLikelihoodExt treeLikelihood) {

        double substitutionRate = 1000.1;//(getHaplotypeCount()*getHaplotypeLength()) ;
        double damageRate = 0;
        SiteModel siteModel = treeLikelihood.getSiteModel();
        SubstitutionModel substitutionModel = siteModel.getSubstitutionModel();
        
        int[] initialSequence = aMap.getConsensusSequenceState();
        StringBuffer buffer = new StringBuffer();
        for (int i = 0; i < initialSequence.length; i++) {
            buffer.append(Nucleotides.INSTANCE.getChar(    initialSequence[i] ));
        }
        treeLikelihood.makeDirty();
        treeLikelihood.getLogLikelihood();
        System.err.println(buffer.toString());



//    	Arrays.fill(initialSequence, 1);
        SeqGenExt seqGen = new SeqGenExt(initialSequence, 
                substitutionRate, substitutionModel, siteModel,
                damageRate);
        
        Tree tree = treeLikelihood.getTreeModel();
        jebl.evolution.alignments.Alignment jeblAlignment = seqGen.simulate(tree);
        List<Sequence> sequenceList = jeblAlignment.getSequenceList();
        for (int j = 0; j < sequenceList.size(); j++) {
			

//        	System.err.println(getHaplotypeString(j));
			System.out.println(sequenceList.get(j).getString());
			Haplotype haplotype = getHaplotype(j);
			haplotype.setSequenceString(sequenceList.get(j).getString());
			System.out.println(haplotype.getSequenceString());
			System.out.println();
		}
//        SimpleAlignment
		
		
	}

	public void simulateSequence(double errorRate, SiteModel siteModel, SubstitutionModel substitutionModel,
			TreeModel treeModel) {

        double substitutionRate = errorRate/(getHaplotypeCount()*getHaplotypeLength()) ;
        System.err.println(substitutionRate);
        double damageRate = 0;
//        SiteModel siteModel = treeLikelihood.getSiteModel();
//        SubstitutionModel substitutionModel = siteModel.getSubstitutionModel();
        
        int[] initialSequence = aMap.getConsensusSequenceState();
        StringBuffer buffer = new StringBuffer();
        for (int i = 0; i < initialSequence.length; i++) {
            buffer.append(Nucleotides.INSTANCE.getChar(    initialSequence[i] ));
        }

        SeqGenExt seqGen = new SeqGenExt(initialSequence, 
                substitutionRate, substitutionModel, siteModel, damageRate);
        
        jebl.evolution.alignments.Alignment jeblAlignment = seqGen.simulate(treeModel);
        List<Sequence> sequenceList = jeblAlignment.getSequenceList();
        for (int j = 0; j < sequenceList.size(); j++) {
			Haplotype haplotype = getHaplotype(j);
			haplotype.setSequenceString(sequenceList.get(j).getString());
		}
		
	}
	
	public static HaplotypeModel factory(Alignment shortReads, Alignment trueAlignment){
		
		HaplotypeModel haplotypeModel = new HaplotypeModel(trueAlignment);
		return haplotypeModel;
	}




	public static HaplotypeModel duplicateHaplotypeModel(
			HaplotypeModel oldModel) {
		HaplotypeModel haplotypeModel = new HaplotypeModel(oldModel.getHaplotypeCount(), oldModel.getHaplotypeLength());

		for (int i = 0; i < oldModel.getHaplotypeCount(); i++) {
			Haplotype haplotype = Haplotype.duplicateHaplotype(oldModel.getHaplotype(i));
			haplotypeModel.setHaplotype(i, haplotype);
		}
		
		return haplotypeModel;

	}	
//XXX not used after this
	public double getLogqFrequency(int oldChar, int newChar){
		return storedLogqMatrix[NUCLEOTIDE_STATES[oldChar]][NUCLEOTIDE_STATES[newChar]];
	}
	
	public double getLogqFrequencyStates(int oldState, int newState){
		return storedLogqMatrix[oldState][newState];
	}
	
	private void checkFrequencyParameter(Parameter frequency) {

		for (int i = 0; i < storedFrequency.length; i++) {
			if(storedFrequency[i]!= frequency.getParameterValue(i)){

				for (int j = i; j < storedFrequency.length; j++) {
					storedFrequency[j] = frequency.getParameterValue(j);
					logFreq[j] = Math.log(storedFrequency[j]);
				}			
				storedCumSumFrequency[0] = storedFrequency[0];
				storedCumSumFrequency[1] = storedCumSumFrequency[0]+storedFrequency[1];
				storedCumSumFrequency[2] = storedCumSumFrequency[1]+storedFrequency[2];
//				storedCumSumFrequency[2] = storedCumSumFrequency[1]+storedFrequency[2];

				for (int j = 0; j < logFreq.length; j++) {
					for (int k = j+1; k < logFreq.length; k++) {
						storedLogqMatrix[j][k] = logFreq[j]-logFreq[k];
						storedLogqMatrix[k][j] = -storedLogqMatrix[j][k];
					}
//					System.out.println(Arrays.toString(storedLogqMatrix[j]));
				}
				
				break;
			}
		}
	}

	private void checkFrequencyParameterTune(Parameter frequency) {

		for (int i = 0; i < storedFrequency.length; i++) {
			if(storedFrequency[i]!= frequency.getParameterValue(i)){
				System.out.println(i);
long time1 = System.currentTimeMillis();
for (int t = 0; t < 1e7; t++) {
	
				for (int j = i; j < storedCumSumFrequency.length; j++) {
					storedFrequency[j] = frequency.getParameterValue(j);
					logFreq[j] = Math.log(storedFrequency[j]);
				}			
//				storedFrequency = frequency.getParameterValues();
				storedCumSumFrequency[0] = storedFrequency[0];
				storedCumSumFrequency[1] = storedCumSumFrequency[0]+storedFrequency[1];
				storedCumSumFrequency[2] = storedCumSumFrequency[1]+storedFrequency[2];
				// Too short for a loop?
				// for (int j = 1; j < INDEX_OF_LAST_VALID_CHARS; j++) {
				// storedCumSumFrequency[j] = storedCumSumFrequency[j - 1]
				// + storedFrequency[j];
				// }
				
//				double[] logFreq = new double[4];
//				for (int j = 0; j < storedCumSumFrequency.length; j++) {
//					logFreq[j] = Math.log(storedFrequency[j]);
//				}
//				System.err.println(Arrays.toString(storedFrequency));
				for (int j = 0; j < storedFrequency.length; j++) {
					for (int k = j+1; k < storedFrequency.length; k++) {
//						storedLogqMatrix[j][k] = Math.log(storedFrequency[j]/storedFrequency[k]);
						storedLogqMatrix[j][k] = logFreq[j]-logFreq[k];
						storedLogqMatrix[k][j] = -storedLogqMatrix[j][k];
					}
//					System.out.println(Arrays.toString(storedLogqMatrix[j]));
				}
//				System.out.println();
//				for (int j = 0; j < storedFrequency.length; j++) {
//					for (int k = 0; k < storedFrequency.length; k++) {
////						storedLogqMatrix[j][k] = Math.log(storedFrequency[j]/storedFrequency[k]);
//						storedLogqMatrix[j][k] = logFreq[j]-logFreq[k];
//					}
////					System.out.println(Arrays.toString(storedLogqMatrix[j]));
//				}

//				System.out.println(i +"\t"+ Arrays.toString(storedFrequency));
//				System.err.println(Arrays.toString(storedCumSumFrequency));
}
long time2 = System.currentTimeMillis();

System.out.println((time2 - time1) + "\t");
				
				break;
			}
		}
	}
	
	private double[] logFreq = new double[4];
	private double[] storedFrequency = new double[4];
	private double[] storedCumSumFrequency = new double[STATE_COUTN];
	private double[][] storedLogqMatrix = new double[4][4];

	@Deprecated AlignmentMapping aMap;
	
}