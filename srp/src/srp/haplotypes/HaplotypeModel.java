package srp.haplotypes;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jebl.evolution.sequences.Sequence;
import srp.dr.ext.SeqGenExt;
import srp.dr.ext.TreeLikelihoodExt;
import srp.evolution.OperationType;
import srp.evolution.haplotypes.old.OldHapOperation;
import srp.evolution.haplotypes.old.OldHapSwapInfo;
import srp.evolution.haplotypes.old.OldHaplotypeModel;
import srp.evolution.shortreads.AlignmentMapping;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;
import dr.util.NumberFormatter;

public class HaplotypeModel extends AbstractHaplotypeModel  {
	
	private static final String MODEL_NAME = "HaplotypeModel";
	private static final long serialVersionUID = -5057514703825711955L;

	private boolean DEBUG = false;
	
//	private static final int HAP_INDEX = OldHapSwapInfo.SWAPBASE_HAP_INDEX;
//	private static final int POS_INDEX = OldHapSwapInfo.SWAPBASE_POS_INDEX;
//	private static final int NEW_CHAR_INDEX = OldHapSwapInfo.SWAPBASE_NEW_CHAR_INDEX;
//	private static final int OLD_CHAR_INDEX = OldHapSwapInfo.SWAPBASE_OLD_CHAR_INDEX;

	private static final int NUCLEOTIDE_STATES[] = Nucleotides.NUCLEOTIDE_STATES;
	private static final char[] VALID_CHARS = initValidChars4();
	private static final int INDEX_OF_LAST_VALID_CHARS = VALID_CHARS.length-1;
	

//	int haplotypesCount;
	AlignmentMapping aMap;

	
//	protected OperationRecord operationRecord;
	
	private boolean isEdit;
	
	
	private HaplotypeModel(AlignmentMapping aMap){
		super(MODEL_NAME);

		this.aMap = aMap;
		haplotypesLength = this.aMap.getLength();
		
		haplotypes = new ArrayList<Haplotype>();


		
//		storedCumSumFrequency[INDEX_OF_LAST_VALID_CHARS]=1;
	}
	private HaplotypeModel(Alignment shortReads) {
		this(new AlignmentMapping(shortReads));
	}


	private static char[] initValidChars4() {
		char[] validChar = new char[4];
		System.arraycopy(Nucleotides.NUCLEOTIDE_CHARS, 0, validChar, 0, validChar.length);
		return validChar;
	}

	private void addHaplotype(Haplotype haplotype) {

	    haplotypes.add(haplotype);
	
	}

	private void initSeqs(int hapCount){
		
		char[] temp = new char[haplotypesLength];
		Arrays.fill(temp, DataType.GAP_CHARACTER);
		String tempSeq = String.valueOf(temp);
		
		for (int i = 0; i < hapCount; i++) {
			Taxon t = new Taxon(TAXON_PREFIX+i); 
			Haplotype haplotype = new Haplotype(t, tempSeq);
			addHaplotype(haplotype);
			randomHaplotype(i);

		}

//
//
//		int i=0;
//		Taxon t = new Taxon(TAXON_PREFIX+i); 
//		Haplotype haplotype = new Haplotype(t, tempSeq);
//		addHaplotype(haplotype);
//		randomHaplotype(i);
//
//		for (i = 1; i < hapCount; i++) {
//			t = new Taxon(TAXON_PREFIX+i); 
//			Haplotype haplotype2 = new Haplotype(t, haplotype.getSequenceString());
//			addHaplotype(haplotype2);
//		}


	}
	

	private void randomHaplotype(int hapIndex) {
		Haplotype haplotype = haplotypes.get(hapIndex);
		for (int i = 0; i < haplotypesLength; i++) {
			//TODO: Which is the "good" starting point
//			char newChar = (char) aMap.nextBaseAt(i);
			char newChar = (char) aMap.nextBaseEqualFreq();
			haplotype.setCharAt(i, newChar);
			
		}
	}
	
	public HaplotypeModel(AlignmentMapping aMap, int hapCount) {
		this(aMap);		
		initSeqs(hapCount);
		
	}

	public HaplotypeModel(AlignmentMapping aMap, Alignment trueAlignment) {
		this(aMap);

		for (int i = 0; i < trueAlignment.getSequenceCount(); i++) {
			Haplotype haplotype = new Haplotype(trueAlignment.getSequence(i));
			addHaplotype(haplotype);
		}
	}
	
	public HaplotypeModel(Alignment shortReads, int hapCount) {
		this(new AlignmentMapping(shortReads), hapCount);

	}	
	
/////////////////////////////////////

	
//	public abstract int getHaplotypeCount();
//	public abstract AbstractSpectrum getHaplotype(int i);
//	public abstract void removeSpectrum(int i);
//	public abstract String getSpectrumString(int i);

	
//	public double[] getSpecturmFrequencies(int spectrumIndex, int i) {
//			return getHaplotype(spectrumIndex).getFrequenciesAt(i);
//	}
//
//	public void resetHaplotypeOperation () {
//		operationRecord.setOperation(OperationType.FULL);
//	}
//
//	public OperationRecord getHaplotypeOperationRecord() {
//		return operationRecord;
//	}
//
//	public OperationType getHaplotypeOperation () {
//		return operationRecord.getOperation();
//	}
//
//	public void setHaplotypeOperationRecord(OperationType op,
//			int[] twoSpectrumIndex, int[] swapPositionIndex) {
//
//		operationRecord.setRecord(op, twoSpectrumIndex,
//				swapPositionIndex);
//	}
//
//	public void setHaplotypeOperationRecord(OperationType op, int siteIndex,
//			double... delta) {
//		operationRecord.setRecord(op, siteIndex, delta);
//	}
//
//	public void setHaplotypeOperationRecord(OperationType op,
//			int spectrumIndex, int siteIndex) {
//		operationRecord.setRecord(op, spectrumIndex, siteIndex);
//	}
//
//	public void setHaplotypeOperationRecord(OperationType op,
//			int spectrumIndex, int siteIndex, double delta) {
//		operationRecord.setRecord(op, spectrumIndex, siteIndex);
//	}
//
//	public void setHaplotypeOperationRecord(OperationType op,
//			int spectrumIndex, int[] siteIndexs, double... delta) {
//		operationRecord.setRecord(op, spectrumIndex, siteIndexs, delta);
//
//	}
//
//	public void setHaplotypeOperationRecord(OperationType op,
//			int[] spectrumIndexs, int siteIndex) {
//		// subcolumn
//		operationRecord.setRecord(op, spectrumIndexs, siteIndex);
//
//	}


	
	
	
	
////////////////////////////////////
	

	
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
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		System.err.println("Call handleModelChangedEvent");
		
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		System.err.println("Call handleVariableChangedEvent");
	}

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
	
	
	public AlignmentMapping getAlignmentMapping() {
		return aMap;
	}

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
	private double[] storedCumSumFrequency = new double[INDEX_OF_LAST_VALID_CHARS];
	private double[][] storedLogqMatrix = new double[4][4];


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
			System.out.println(getHaplotypeString(j));
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
////////////////////////////////////////////
/// remove later, will be replaced in operator 
	@Deprecated private OldHapSwapInfo swapInfo = new OldHapSwapInfo();
	@Deprecated int[] swapBaseRecord = new int[4];
	@Deprecated int HAP_INDEX = 0;
	@Deprecated int POS_INDEX = 1;
	@Deprecated int NEW_CHAR_INDEX = 2;
	@Deprecated int OLD_CHAR_INDEX = 3;
	
	@Deprecated
	public void swapHaplotypeColumn(int[] posChar){

		int[] allOldChars = new int[getHaplotypeCount()];
			
		for (int i = 0; i < getHaplotypeCount(); i++) {
			Haplotype haplotype = haplotypes.get(i);
			allOldChars[i] = haplotype.replaceCharAt(posChar[0], posChar[1]);
		}
			
		storeOperationRecord(OldHapOperation.SWAPCOLUMN, posChar, allOldChars);


	}
	@Deprecated
	public int[] swapHaplotypeSingleBase(OldHapOperation op, int[] posChar){

		swapBaseRecord[POS_INDEX] = posChar[0];
		swapBaseRecord[NEW_CHAR_INDEX] = posChar[1];

		swapBaseRecord[HAP_INDEX] = MathUtils.nextInt(getHaplotypeCount());

		Haplotype haplotype = haplotypes.get(swapBaseRecord[HAP_INDEX]);
		swapBaseRecord[OLD_CHAR_INDEX] = haplotype.replaceCharAt(
				swapBaseRecord[POS_INDEX], swapBaseRecord[NEW_CHAR_INDEX]);

		storeOperationRecord(op, swapBaseRecord);
		return swapBaseRecord;

	}

	@Deprecated
	public void swapHaplotypeMultiBases(OldHapOperation op, int[] allNewChars){
		int hapIndex = MathUtils.nextInt(getHaplotypeCount());
		
		Haplotype haplotype = haplotypes.get(hapIndex);
		
		int[] allOldChars = new int[allNewChars.length];
		for (int i = 0; i < haplotypesLength; i++) {
			int newChar = allNewChars[i];
			if(newChar>0){
				allOldChars[i] = haplotype.replaceCharAt(i, newChar);
			}
		}
		storeOperationRecord(op, new int[]{hapIndex}, allNewChars, allOldChars);
//		swapInfo.storeOperation(op, new int[]{hapIndex}, allPosChars);
		
	}
	@Deprecated
	public void swapHaplotypeMultiBases(OldHapOperation op, int[][] allPosChars){
		int hapIndex = MathUtils.nextInt(getHaplotypeCount());
		
		Haplotype haplotype = haplotypes.get(hapIndex);

		for (int i = 0; i < haplotypesLength; i++) {
			int newChar = allPosChars[0][i];
			if(newChar>0){
				allPosChars[1][i] = haplotype.replaceCharAt(i, newChar);
			}
		}
		storeOperationRecord(op, hapIndex, allPosChars);
//		storeOperationRecord(op, new int[] {hapIndex}, allPosChars[0], allPosChars[1]);
//		swapInfo.storeOperation(op, new int[]{hapIndex}, allPosChars);
		
	}
	@Deprecated
	private char setNewCharFromFrequency(){
		
		double d = MathUtils.nextDouble();
		
		for (int i = 0; i < INDEX_OF_LAST_VALID_CHARS; i++) {
			if (d <= storedCumSumFrequency[i]) {
				return VALID_CHARS[i];
			}
		}
		return VALID_CHARS[INDEX_OF_LAST_VALID_CHARS];
	}
	@Deprecated
	public double swapNextDiffBaseFrequency(OldHapOperation op, Parameter frequency) {
	
		checkFrequencyParameter(frequency);
		
//		int[] swapBaseRecord = new int[4];
		swapBaseRecord[HAP_INDEX] = MathUtils.nextInt(getHaplotypeCount());
		swapBaseRecord[POS_INDEX] = MathUtils.nextInt(getHaplotypeLength());

		Haplotype haplotype = haplotypes.get(swapBaseRecord[HAP_INDEX]);
		swapBaseRecord[OLD_CHAR_INDEX] = haplotype.getChar(swapBaseRecord[POS_INDEX]);
		
		do{
			swapBaseRecord[NEW_CHAR_INDEX] = setNewCharFromFrequency();
		}while(swapBaseRecord[OLD_CHAR_INDEX]==swapBaseRecord[NEW_CHAR_INDEX]);

		haplotype.setCharAt(swapBaseRecord[POS_INDEX], (char) swapBaseRecord[NEW_CHAR_INDEX]);
		
		storeOperationRecord(op, swapBaseRecord);
//			return swapBaseRecord;
//			swapInfo.storeOperation(op, swapInfoArray);
//			int newChar = NUCLEOTIDE_STATES[swapRecord[NEW_CHAR_INDEX]];
//			int oldChar = NUCLEOTIDE_STATES[swapRecord[OLD_CHAR_INDEX]];
		double logq = getLogqFrequency(swapBaseRecord[OLD_CHAR_INDEX], swapBaseRecord[NEW_CHAR_INDEX]);
		
		return logq;
	}
	@Deprecated
	public int[] getNextBaseFrequency(Parameter frequency) {
	
		checkFrequencyParameter(frequency);
	
		int[] tempPosChar = new int[2];
		tempPosChar[0] = MathUtils.nextInt(getHaplotypeLength());
	
		double d = MathUtils.nextDouble();
	
//		for (int i = 0; i < INDEX_OF_LAST_VALID_CHARS; i++) {
//			if (d <= storedCumSumFrequency[i]) {
//				tempPosChar[1] = VALID_CHARS[i];
//				return tempPosChar;
//			}
//		}
		tempPosChar[1] = setNewCharFromFrequency();
		return tempPosChar;
	}

	//	private int replaceHaplotypeCharAt(int hapIndex, int pos, int newChar){
//		
//		int oldChar = haplotypes.get(hapIndex).replaceCharAt(pos, (char) newChar);
//		return oldChar;
//		
//	}
//	private static int replaceHaplotypeCharAt(Haplotype haplotype, int pos, int newChar){
//		int oldChar = haplotype.replaceCharAt(pos, (char) newChar);
//		return oldChar;
//		
//	}
	@Deprecated
	private void resetHaplotypeToOldChar(int[] swapRecord){
		Haplotype haplotype = haplotypes.get(swapRecord[0]);
		haplotype.replaceCharAt(swapRecord[1], swapRecord[3]);

//		swapHaplotypeCharAt(swapArray[0], swapArray[1], swapArray[3]);
	}
	
//	private int swapHaplotypeCharAt(int hapIndex, int[] posChar){
//		int oldChar = getHaplotype(hapIndex).getChar(posChar[0]); 
//		getHaplotype(hapIndex).setCharAt(posChar[0], (char) posChar[1]);
//		return oldChar;
//		
//	}
//	
//	
	@Deprecated
	public OldHapSwapInfo getSwapInfo() {
		return swapInfo;
	}
	@Deprecated
	public void storeOperationRecord(OldHapOperation op, int[]... opRecord){
		swapInfo.storeOperation(op, opRecord);
	}
	@Deprecated
	private void storeOperationRecord(OldHapOperation op, int hapIndex,
			int[][] opRecord) {
		swapInfo.storeOperation(op, opRecord[0], opRecord[1]);
		swapInfo.storeHapIndex(hapIndex);
	}


	public static HaplotypeModel factory(Alignment shortReads, Alignment trueAlignment){
		
		AlignmentMapping alignmentMapping = new AlignmentMapping(shortReads);
		HaplotypeModel haplotypeModel = new HaplotypeModel(alignmentMapping, trueAlignment);
		return haplotypeModel;
	}	
}