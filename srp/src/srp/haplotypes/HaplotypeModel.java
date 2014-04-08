package srp.haplotypes;


import java.util.List;

import jebl.evolution.sequences.Sequence;
import srp.dr.ext.SeqGenExt;
import srp.dr.ext.TreeLikelihoodExt;
import srp.evolution.OperationType;
import srp.evolution.shortreads.AlignmentMapping;
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
	
	private static final String MODEL_NAME = "HaplotypeModel";
	private static final long serialVersionUID = -5057514703825711955L;

	private boolean DEBUG = false;

	private static final int NUCLEOTIDE_STATES[] = Nucleotides.NUCLEOTIDE_STATES;
//	private static final char[] VALID_CHARS = initValidChars4();
	private static final int STATE_COUTN = DATA_TYPE.getStateCount();
	


	private void initHaplotypes() {
		for (int i = 0; i < haplotypeCount; i++) {
			Taxon taxon = new Taxon(TAXON_PREFIX + i);
			Haplotype haplotype = new Haplotype(taxon, haplotypeLength);
			setHaplotype(i, haplotype);
		}

	}
	

	
	public HaplotypeModel(int hapCount, int hapLength) {
		super(MODEL_NAME, hapCount, hapLength);
		initHaplotypes();
//		this.haplotypeCount = hapCount;
//		this.haplotypeLength = hapLength;

//		haplotypes = new ArrayList<Haplotype>();

		
	}

	public HaplotypeModel(Alignment trueAlignment) {
		this(trueAlignment.getSequenceCount(), trueAlignment.getSiteCount());

		for (int i = 0; i < trueAlignment.getSequenceCount(); i++) {
			Haplotype haplotype = new Haplotype(trueAlignment.getSequence(i));
			setHaplotype(i, haplotype);
		}
	}
	
//	public HaplotypeModel(Alignment shortReads, int hapCount) {
//		this(new AlignmentMapping(shortReads), hapCount);
//
//	}	
	

	
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



    
	
//	@Override
//	public void fireModelChanged(){
////		for (TreeChangedEvent treeChangedEvent : treeChangedEvents) {
//		listenerHelper.fireModelChanged(this);//, treeChangedEvent);
////		}
////		treeChangedEvents.clear();
//	}


	@Override
	protected void acceptState() {
		//Do nothing
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
			if(DEBUG){
				System.out.println("StoreState in SpectrumAlignment:\t"+operation);
			}
			break;
		case FULL:
			if(DEBUG){
				System.out.println("StoreState in SpectrumAlignment:\t"+operation);
			}
			for (int i = 0; i < getHaplotypeCount(); i++) {
				haplotype = getHaplotype(i);
				for (int s = 0; s < getHaplotypeLength(); s++) {
					haplotype.storeState(s);
				}
			}
			break;
//		case SINGLE_DELTA:
//			spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
//			siteIndex = spectrumOperationRecord.getAllSiteIndexs()[0];
//
//			spectrum = getHaplotype(spectrumIndex);
//			spectrum.setStoreSiteIndex(siteIndex);
//			spectrum.storeState();
//			break;
		case COLUMN:

//			spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			siteIndex = operationRecord.getSingleIndex();
			for (int i = 0; i < getHaplotypeCount(); i++) {
				haplotype = getHaplotype(i);
				haplotype.storeState(siteIndex);
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
		case RECOMBINATION:

//			System.err.println("store alignment recombination");
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
			throw new IllegalArgumentException("Unknown operation type: "+operation +"\tin"+SpectrumAlignmentModel.class.getSimpleName() );
			
		}

//		 
//		long time2 = System.currentTimeMillis();
//		time += (time2-time1);
		 
	}
	@Override
	protected void restoreState() {
//		long time1 = System.currentTimeMillis();
		
		OperationType operation = operationRecord.getOperation();
		int spectrumIndex;
		int siteIndex;
		Haplotype haplotype;
		int[] siteIndexs;
//			System.out.println(op);
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
				for (int s = 0; s < getHaplotypeLength(); s++) {
					haplotype.restoreState(s);
				}
			}
			break;
//			
//				spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
//				siteIndex = spectrumOperationRecord.getAllSiteIndexs()[0];
////				System.err.println(spectrumIndex +"\t"+ siteIndex +"\t"+ Arrays.toString(getSpecturmFrequencies(spectrumIndex,
////						siteIndex)));
//				spectrum = getHaplotype(spectrumIndex);
//				spectrum.setStoreSiteIndex(siteIndex);
//				spectrum.restoreState();
////				System.err.println("after restore\t"+spectrumIndex +"\t"+ siteIndex +"\t"+ Arrays.toString(getSpecturmFrequencies(spectrumIndex,
////						siteIndex)));
//				break;
		case COLUMN:

//				spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			siteIndex = operationRecord.getSingleIndex();
//				System.err.println(spectrumIndex +"\t"+ siteIndex +"\t"+ Arrays.toString(getSpecturmFrequencies(spectrumIndex,
//						siteIndex)));
			for (int i = 0; i < getHaplotypeCount(); i++) {
				haplotype = getHaplotype(i);
				haplotype.restoreState(siteIndex);
			}
			break;
		case MULTI:

			spectrumIndex = operationRecord.getSpectrumIndex();
			siteIndexs = operationRecord.getAllSiteIndexs();
			haplotype = getHaplotype(spectrumIndex);
			for (int s = 0; s < siteIndexs.length; s++) {
				haplotype.restoreState(s);
			}
			break;
		case SINGLE:

			spectrumIndex = operationRecord.getSpectrumIndex();
			siteIndex = operationRecord.getSingleIndex();
			haplotype = getHaplotype(spectrumIndex);

			haplotype.restoreState(siteIndex);

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
//		long time2 = System.currentTimeMillis();
//		time += (time2-time1);
	}
	
	
//	public AlignmentMapping getAlignmentMapping() {
//		return aMap;
//	}

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
}