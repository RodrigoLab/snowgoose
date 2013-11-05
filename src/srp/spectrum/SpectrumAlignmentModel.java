package srp.spectrum;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import srp.haplotypes.AbstractHaplotypeModel;
import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.Haplotype;
import srp.haplotypes.Operation;
import srp.haplotypes.SwapInfo;

import jebl.evolution.sequences.Sequence;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.ext.SeqGenExt;
import dr.ext.TreeLikelihoodExt;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;
import dr.util.NumberFormatter;



public class SpectrumAlignmentModel extends AbstractSpectrumAlignmentModel  {
	
	private static final long serialVersionUID = 3458306765918357829L;
	private static final String MODEL_NAME = "SpectrumModel";

	
	private static final int HAP_INDEX = SwapInfo.SWAPBASE_HAP_INDEX;
	private static final int POS_INDEX = SwapInfo.SWAPBASE_POS_INDEX;
	private static final int NEW_CHAR_INDEX = SwapInfo.SWAPBASE_NEW_CHAR_INDEX;
	private static final int OLD_CHAR_INDEX = SwapInfo.SWAPBASE_OLD_CHAR_INDEX;

	private static final int NUCLEOTIDE_STATES[] = Nucleotides.NUCLEOTIDE_STATES;
	private static final char[] VALID_CHARS = initValidChars4();
	private static final int INDEX_OF_LAST_VALID_CHARS = VALID_CHARS.length-1;
	
//	private static final DataType DATA_TYPE = Nucleotides.INSTANCE;
	public static final char GAP = '-';
	public static final String TAXON_PREFIX = "taxa_";

//	int haplotypesCount;
	AlignmentMapping aMap;
	private SwapInfo swapInfo = new SwapInfo();
	
	private boolean isEdit;
	
	int[] swapBaseRecord = new int[4];
	
	private SpectrumAlignmentModel(AlignmentMapping aMap){
		super(MODEL_NAME);

		this.aMap = aMap;
		spectrumLength = this.aMap.getLength();
		
		spectrumList = new ArrayList<Spectrum>();
		setDataType(Nucleotides.INSTANCE);

		
//		storedCumSumFrequency[INDEX_OF_LAST_VALID_CHARS]=1;
	}
	private SpectrumAlignmentModel(Alignment shortReads) {
		this(new AlignmentMapping(shortReads));
	}


	private static char[] initValidChars4() {
		char[] validChar = new char[4];
		System.arraycopy(Nucleotides.NUCLEOTIDE_CHARS, 0, validChar, 0, validChar.length);
		return validChar;
	}

	private void addSpectrum(Spectrum spectrum) {
		spectrum.setDataType(getDataType());
	    spectrumList.add(spectrum);
	
	}

	private void initSeqs(int hapCount){
		
		char[] temp = new char[spectrumLength];
		Arrays.fill(temp, DataType.GAP_CHARACTER);
		String tempSeq = String.valueOf(temp);
		
		for (int i = 0; i < hapCount; i++) {
			Taxon t = new Taxon(TAXON_PREFIX+i); 
			Spectrum spectrum = new Spectrum(spectrumLength);
			spectrum.setTaxon(t);
			addSpectrum(spectrum);
//			randomHaplotype(i);

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
//	
//
//	private void randomHaplotype(int hapIndex) {
//		Spectrum spectrum = spectrumList.get(hapIndex);
//		for (int i = 0; i < spectrumLength; i++) {
//			//TODO: Which is the "good" starting point
////			char newChar = (char) aMap.nextBaseAt(i);
//			char newChar = (char) aMap.nextBaseEqualFreq();
//			haplotype.setCharAt(i, newChar);
//			
//		}
//	}
//	
	public SpectrumAlignmentModel(AlignmentMapping aMap, int hapCount) {
		this(aMap);		
		initSeqs(hapCount);
		
	}

	public SpectrumAlignmentModel(AlignmentMapping aMap, Alignment trueAlignment) {
		this(aMap);

		for (int i = 0; i < trueAlignment.getSequenceCount(); i++) {
			Spectrum spectrum = new Spectrum(trueAlignment.getSequence(i));
			addSpectrum(spectrum);
		}
	}
	
	public SpectrumAlignmentModel(Alignment shortReads, int hapCount) {
		this(new AlignmentMapping(shortReads), hapCount);

	}	
	

	

	private char setNewCharFromFrequency(){
		
		double d = MathUtils.nextDouble();
		
		for (int i = 0; i < INDEX_OF_LAST_VALID_CHARS; i++) {
			if (d <= storedCumSumFrequency[i]) {
				return VALID_CHARS[i];
			}
		}
		return VALID_CHARS[INDEX_OF_LAST_VALID_CHARS];
	}
	

	public int[] getNextBaseFrequency(Parameter frequency) {
	
		checkFrequencyParameter(frequency);
	
		int[] tempPosChar = new int[2];
		tempPosChar[0] = MathUtils.nextInt(getSpectrumLength());
	
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
	private void resetHaplotypeToOldChar(int[] swapRecord){
//		Haplotype haplotype = spectrumList.get(swapRecord[0]);
//		haplotype.replaceCharAt(swapRecord[1], swapRecord[3]);

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

	public SwapInfo getSwapInfo() {
		return swapInfo;
	}
//	public Operation getOperation() {
//		return swapInfo.getOperation();
//	}
	public void storeOperationRecord(Operation op, int[]... opRecord){
		swapInfo.storeOperation(op, opRecord);
	}
	
	private void storeOperationRecord(Operation op, int hapIndex,
			int[][] opRecord) {
		swapInfo.storeOperation(op, opRecord[0], opRecord[1]);
		swapInfo.storeHapIndex(hapIndex);
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
        for (int i = 0; i < getSpectrumCount(); i++) {
            String name = formatter.formatToFieldWidth(getTaxon(i).getId(), 10);
            buffer.append(">" + name + "\n");
            buffer.append(getAlignedSequenceString(i) + "\n");
        }

        return buffer.toString();
    
		
	}



    
	
	@Override
		public void fireModelChanged(){
	//		for (TreeChangedEvent treeChangedEvent : treeChangedEvents) {
			listenerHelper.fireModelChanged(this);//, treeChangedEvent);
	//		}
	//		treeChangedEvents.clear();
		}

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

		// System.err.println("Call storeState");
		swapInfo.storeOperation(Operation.NONE);
	}

	@Override
	protected void restoreState() {
		// System.err.println("Call restoreState - haplotypeModel");
		reject();
		swapInfo.storeOperation(Operation.NONE);
	}

	//	private int swapHaplotypeCharAt(int hapIndex, int[] posChar){
	//		int oldChar = getHaplotype(hapIndex).getChar(posChar[0]); 
	//		getHaplotype(hapIndex).setCharAt(posChar[0], (char) posChar[1]);
	//		return oldChar;
	//		
	//	}
	//	
	//	

	public void reject() {
		//TODO redo
		Operation op = swapInfo.getOperation();
//		int[] temp;
		switch (op) {
		
		case NONE:
			break;
		case PASS:
			break;
//		case SWAPSINGLE:
//
//			int[] temp = swapInfo.getSwapInfoSWAPBASE();
//			//			long time1 = System.currentTimeMillis();
////			for (int i = 0; i < 1e9; i++) {
//			resetHaplotypeToOldChar(temp);
////		replaceHaplotypeCharAt(haplotypes.get(temp[0]), temp[2], temp[3]);
////			}
////			long time2 = System.currentTimeMillis();
////
////			System.out.println("Single: "+(time2 - time1) + "\t");
////			System.exit(0);
//			break;

//		case SWAPMULTI:
//			
//			int hapIndex = swapInfo.getHapIndex();
//			int[][] swapMulti = swapInfo.getSwapInfoSWAPMULTI();
////			int[] allNewChars = swapMulti[0];
//			int[] allOldChars = swapMulti[1];
//			
//			
//			Haplotype haplotype = spectrumList.get(hapIndex);
//			for (int i = 0; i < allOldChars.length; i++) {
//				int oldChar = allOldChars[i];
//				if(oldChar>0){
//					haplotype.replaceCharAt(i, oldChar);
//				}
//			}
//			break;
//		
//		case SWAPSECTION:
//			int[] swapSection = swapInfo.getSwapInfoSWAPSECTION();
//			
////			time1 = System.currentTimeMillis();
//////			for (int i = 0; i < 1e8; i++) {
////				haplotypes.get(temp[0]).restoreState();
//////			}
////			time2 = System.currentTimeMillis();
//
////			System.out.println("Multi: "+(time2 - time1) + "\t");
//
//			
//			
//			for (int i = 0; i < 2; i++) {
//				spectrumList.get(swapSection[i]).restoreState();
//			}
//			break;
//		case SWAPCOLUMN:
//			int[][] swapColumn = swapInfo.getSwapInfoSWAPCOLUMN();
//			int pos = swapColumn[0][0];
//			int[] allOldChars2 = swapColumn[1];
//
//			for (int i = 0; i < getHaplotypeCount(); i++) {
//				haplotype = spectrumList.get(i);
//				haplotype.replaceCharAt(pos, allOldChars2[i]);
//			}
//			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: " + op);

		}


	}

	public void startSpectrumOperation(){
		isEdit = true;
	}

	public void endSpectrumOperation(){
		isEdit = false;
		fireModelChanged();
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

	private double[] logFreq = new double[4];
	private double[] storedFrequency = new double[4];
	private double[] storedCumSumFrequency = new double[INDEX_OF_LAST_VALID_CHARS];
	private double[][] storedLogqMatrix = new double[4][4];

	public double[] getSpecturmFrequencies(int sequenceIndex, int i) {
//		Spectrum spectrum = getSpectrum(sequenceIndex);
		
		return getSpectrum(sequenceIndex).getFrequencies(i);
	}
	public void setSpectrum(int index, Spectrum spectrum) {
		// TODO Auto-generated method stub
		spectrumList.set(index, spectrum);
	}


	
	
}