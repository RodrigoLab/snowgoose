package srp.spectrum;


import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.StringTokenizer;

import srp.haplotypes.AlignmentMapping;
import srp.haplotypes.Operation;
import srp.haplotypes.SwapInfo;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
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
	private AlignmentMapping aMap;
	private SwapInfo swapInfo = new SwapInfo();

	private boolean isEdit;
	
	private int[] swapBaseRecord = new int[4];
	private SpectrumOperationRecord spectrumOperationRecord;
	

	
	private SpectrumAlignmentModel(AlignmentMapping aMap){
		super(MODEL_NAME);
		this.aMap = aMap;
		spectrumLength = this.aMap.getLength();
		
		spectrumList = new ArrayList<Spectrum>();
		storedSpectrumList = new ArrayList<Spectrum>();
		setDataType(Nucleotides.INSTANCE);
		spectrumOperationRecord = new SpectrumOperationRecord();
		
//		
//		spectrumLength = this.aMap.getLength();
//		
//		spectrumList = new ArrayList<Spectrum>();
//		storedSpectrumList = new ArrayList<Spectrum>();
//		setDataType(Nucleotides.INSTANCE);
//		spectrumOperationRecord = new SpectrumOperationRecord();
		
//		storedCumSumFrequency[INDEX_OF_LAST_VALID_CHARS]=1;
	}
//	private SpectrumAlignmentModel(Alignment shortReads) {
//		this(new AlignmentMapping(shortReads));
//	}


	private static char[] initValidChars4() {
		char[] validChar = new char[4];
		System.arraycopy(Nucleotides.NUCLEOTIDE_CHARS, 0, validChar, 0, validChar.length);
		return validChar;
	}

	public void addSpectrum(Spectrum spectrum) {
		spectrum.setDataType(getDataType());
		for (int i = 0; i < spectrumLength; i++) {
			spectrum.setStoreSiteIndex(i);
			spectrum.storeState();
		}
//		spectrum.storeState()
	    spectrumList.add(spectrum);
//	    storedSpectrumList.add(new Spectrum(spectrum));
	    
	    
	    
//	    System.err.println(spectrumLength +"\t"+ getSpectrumCount());
//	    Spectrum spectrum2 = spectrumList.get(getSpectrumCount()-1);
//	    for (int i = 0; i < 10; i++) {
//	    	SpectraParameter spectra = spectrum2.getSpectra(i);
//	    	System.out.println(Arrays.toString(spectra.getFrequencies()));
////	    	for (int j = 0; j < 4; j++) {
////				System.out.print(spectra.getf);
////			}
////			System.out.println(spectra );
//		}
//	    System.out.println();
//	    SpectraParameter spectra = spectrum2.getSpectra(1199);
//    	System.out.println(Arrays.toString(spectra.getFrequencies()));
	
	}

	
	public SpectrumAlignmentModel(Alignment shortReads, int hapCount) {
		this(new AlignmentMapping(shortReads), hapCount);
	}
	
	public SpectrumAlignmentModel(AlignmentMapping aMap, int hapCount) {
		this(aMap, hapCount, false);		
	}
	
	public SpectrumAlignmentModel(AlignmentMapping aMap, int hapCount, boolean random) {
		this(aMap);		

		for (int i = 0; i < hapCount; i++) {
			Taxon t = new Taxon(TAXON_PREFIX+i); 
			Spectrum spectrum = new Spectrum(spectrumLength, random);
			spectrum.setTaxon(t);
			addSpectrum(spectrum);
//			randomHaplotype(i);

		}
	}


	
	public SpectrumAlignmentModel(AlignmentMapping aMap, Alignment trueAlignment) {
		this(aMap);

		for (int i = 0; i < trueAlignment.getSequenceCount(); i++) {
			Spectrum spectrum = new Spectrum(trueAlignment.getSequence(i));
			Taxon t = new Taxon(TAXON_PREFIX+i); 
//			Spectrum spectrum = new Spectrum(spectrumLength);
			spectrum.setTaxon(t);
			addSpectrum(spectrum);
		}
	}
	
	
	
	public static SpectrumAlignmentModel duplicateSpectrumAlignmentModel(SpectrumAlignmentModel oldModel){
		
		SpectrumAlignmentModel newSpectrumModel = new SpectrumAlignmentModel(oldModel.getAlignmentMapping());

		for (int i = 0; i < oldModel.getSpectrumCount(); i++) {
			Spectrum spectrum = Spectrum.duplicateSpectrum(oldModel.getSpectrum(i));
			newSpectrumModel.addSpectrum(spectrum);
		}
		
		
		return newSpectrumModel;
	}
	
	public AlignmentMapping getAlignmentMapping() {
		return aMap;	
	}

	public void setSpectrumOperationRecord(SpectrumOperation op,
			int[] twoSpectrumIndex, int[] swapPositionIndex){
		
		spectrumOperationRecord.setRecord(op, twoSpectrumIndex, swapPositionIndex);
	}

	//	public void setSpectrumOperationRecord(SpectrumOperation op, int spectrumIndex,
//			int siteIndex, double... delta){
//
//		spectrumOperationRecord.setRecord(op, spectrumIndex, siteIndex, delta);
//	}

	public void setSpectrumOperationRecord(SpectrumOperation op, int siteIndex,
			double... delta) {
		spectrumOperationRecord.setRecord(op, siteIndex, delta);
	}

	public void setSpectrumOperationRecord(SpectrumOperation op,
			int spectrumIndex, int[] siteIndexs, double... delta) {
		spectrumOperationRecord.setRecord(op, spectrumIndex, siteIndexs, delta);

	}

	public SpectrumOperationRecord getSpectrumOperationRecord() {
		return spectrumOperationRecord;
	}

	public double[] getSpecturmFrequencies(int spectrumIndex, int i) {
	//		Spectrum spectrum = getSpectrum(sequenceIndex);
			
			return getSpectrum(spectrumIndex).getFrequenciesAt(i);
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
		System.err.println("Call handleModelChangedEvent in SpectrumAlignmentModel");
		
	}
	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		System.err.println("Call handleVariableChangedEvent in SpectrumAlignmentModel");
	}
	@Override
	protected void acceptState() {
		//Do nothing
	}
	@Override
	protected void storeState() {

		SpectrumOperation operation = spectrumOperationRecord.getOperation();
		int spectrumIndex;
		int siteIndex;
		int[] siteIndexs;
		Spectrum spectrum;
		switch (operation) {
		case NONE:
//			System.err.println("SpectrumAlignmentModel NONE operation");
			break;
//		case SINGLE_DELTA:
//			spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
//			siteIndex = spectrumOperationRecord.getAllSiteIndexs()[0];
//
//			spectrum = getSpectrum(spectrumIndex);
//			spectrum.setStoreSiteIndex(siteIndex);
//			spectrum.storeState();
//			break;
		case COLUMN_DELTA:
//			spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			siteIndex = spectrumOperationRecord.getColumnIndex();
			for (int i = 0; i < getSpectrumCount(); i++) {
				spectrum = getSpectrum(i);
				spectrum.setStoreSiteIndex(siteIndex);
				spectrum.storeState();
			}
			break;
		case SINGLE_DELTA:
		case MULTI_DELTA:
			spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			siteIndexs = spectrumOperationRecord.getAllSiteIndexs();
			spectrum = getSpectrum(spectrumIndex);
			for (int i = 0; i < siteIndexs.length; i++) {
				spectrum.setStoreSiteIndex(siteIndexs[i]);
				spectrum.storeState();
			}
			break;
		case RECOMBINATION:
//			System.err.println("store alignment recombination");
			int[] twoSpectrums = spectrumOperationRecord.getRecombinationSpectrumIndex();
			int[] twoPositions = spectrumOperationRecord.getRecombinationPositionIndex();
			for (int i : twoSpectrums) {
				spectrum = getSpectrum(i);
//			TODO wrong code below, replace with above, unittest doesn't catch this
//			for (int i = 0; i < twoSpectrums.length; i++) {
//				spectrum = getSpectrum(i);
				for (int s = twoPositions[0]; s < twoPositions[1]; s++) {
					spectrum.setStoreSiteIndex(s);
					spectrum.storeState();
				}
			}
			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: "+operation +"\tin"+SpectrumAlignmentModel.class.getSimpleName() );
			
		}

		 
		 
//		swapInfo.storeOperation(Operation.NONE);
	}
	@Override
	protected void restoreState() {

		SpectrumOperation op = spectrumOperationRecord.getOperation();
		int spectrumIndex;
		int siteIndex;
		Spectrum spectrum;
		int[] siteIndexs;
//			System.out.println(op);
		switch (op) {
		
		case NONE:
			break;
		case PASS:
			break;
//			
//				spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
//				siteIndex = spectrumOperationRecord.getAllSiteIndexs()[0];
////				System.err.println(spectrumIndex +"\t"+ siteIndex +"\t"+ Arrays.toString(getSpecturmFrequencies(spectrumIndex,
////						siteIndex)));
//				spectrum = getSpectrum(spectrumIndex);
//				spectrum.setStoreSiteIndex(siteIndex);
//				spectrum.restoreState();
////				System.err.println("after restore\t"+spectrumIndex +"\t"+ siteIndex +"\t"+ Arrays.toString(getSpecturmFrequencies(spectrumIndex,
////						siteIndex)));
//				break;
		case COLUMN_DELTA:
//				spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			siteIndex = spectrumOperationRecord.getColumnIndex();
//				System.err.println(spectrumIndex +"\t"+ siteIndex +"\t"+ Arrays.toString(getSpecturmFrequencies(spectrumIndex,
//						siteIndex)));
			for (int i = 0; i < getSpectrumCount(); i++) {
				spectrum = getSpectrum(i);
				spectrum.setStoreSiteIndex(siteIndex);
				spectrum.restoreState();
			}
			break;
		case SINGLE_DELTA:
		case MULTI_DELTA:
			spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			siteIndexs = spectrumOperationRecord.getAllSiteIndexs();
			spectrum = getSpectrum(spectrumIndex);
			for (int i = 0; i < siteIndexs.length; i++) {
				spectrum.setStoreSiteIndex(siteIndexs[i]);
				spectrum.restoreState();
			}
			break;
			
		case RECOMBINATION:
//			System.err.println("Restore alignment recombination");
			int[] twoSpectrums = spectrumOperationRecord.getRecombinationSpectrumIndex();
			int[] twoPositions = spectrumOperationRecord.getRecombinationPositionIndex();
			for (int i : twoSpectrums) {
				spectrum = getSpectrum(i);
//			TODO wrong code below, replace with above, unittest doesn't catch this
//			for (int i = 0; i < twoSpectrums.length; i++) {
//				spectrum = getSpectrum(i);
				for (int s = twoPositions[0]; s < twoPositions[1]; s++) {
					spectrum.setStoreSiteIndex(s);
					spectrum.restoreState();
				}
			}
			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: " + op);

		}

	
	}
	public void startSpectrumOperation(){
//		System.err.println("\n!!!startSpectrumOperation");
		isEdit = true;
	}
	public void endSpectrumOperation(){
		isEdit = false;
//		System.err.println("!!!!!!!endSpectrumOperation, fireModelCHanged");
		fireModelChanged();
	}
//////////////////////////////////////////////////////////////////////////
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
	@Deprecated
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
	@Deprecated
	public SwapInfo getSwapInfo() {
		return swapInfo;
	}
//	public Operation getOperation() {
//		return swapInfo.getOperation();
//	}
	@Deprecated
	public void storeOperationRecord(Operation op, int[]... opRecord){
		swapInfo.storeOperation(op, opRecord);
	}
	@Deprecated
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
//            buffer.append(getAlignedSequenceString(i) + "\n");
        }

        return buffer.toString();
    
		
	}



    
	
	@Deprecated
	public double getLogqFrequency(int oldChar, int newChar){
		return storedLogqMatrix[NUCLEOTIDE_STATES[oldChar]][NUCLEOTIDE_STATES[newChar]];
	}
	@Deprecated
	public double getLogqFrequencyStates(int oldState, int newState){
		return storedLogqMatrix[oldState][newState];
	}
	@Deprecated
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


	public String diagnostic(){

//		int spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
//		int siteIndex = spectrumOperationRecord.getAllSiteIndexs()[0];
//
//		Spectrum spectrum = getSpectrum(spectrumIndex);
//		spectrum.setStoreSiteIndex(siteIndex);
//		SpectraParameter spectra = spectrum.getSpectra(siteIndex);
//		return siteIndex +"\t"+ spectrumIndex +"\t"+ spectra.diagnostic();
//		.restoreState();
		
		int[] twoSpectrums = spectrumOperationRecord.getRecombinationSpectrumIndex();
		int[] twoPositions = spectrumOperationRecord.getRecombinationPositionIndex();
		for (int i = 0; i < twoSpectrums.length; i++) {
		
		}
		return Arrays.toString(twoSpectrums) +"\t"+ Arrays.toString(twoPositions);
//		return "";
	}
	
	public static void compareTwoSpectrumModel(
			SpectrumAlignmentModel spectrumModel,
			SpectrumAlignmentModel spectrumModel2) {


		for (int i = 0; i < spectrumModel.getSpectrumCount(); i++) {
			Spectrum spectrum = spectrumModel.getSpectrum(i);
			Spectrum newSpectrum = spectrumModel2.getSpectrum(i);
			for (int j = 0; j < spectrum.getLength(); j++) {
				double[] frequencies = spectrum.getFrequenciesAt(j);
				for (int k = 0; k < frequencies.length; k++) {
					if(frequencies[k] != newSpectrum.getFrequency(j, k)){
						System.err.println("DIFFMODEL"+i +"\t"+ j +"\t"+ k +"\t"+ Arrays.toString(frequencies) +"\t"+ Arrays.toString(newSpectrum.getFrequenciesAt(j)));
					}
				}
			}
			
		}
		SpectrumOperationRecord record = spectrumModel.getSpectrumOperationRecord();
		int spec = record.getSpectrumIndex();
		int site = record.getAllSiteIndexs()[0];
		System.out.println("compare two models");
		System.out.println(Arrays.toString(spectrumModel.getSpectrum(spec).getFrequenciesAt(
				site)));
		System.out.println(Arrays.toString(spectrumModel2.getSpectrum(spec).getFrequenciesAt(
				site)));
		
	}
	public void removeSpectrum(int i) {
		spectrumList.remove(i);
		
	}
	public StringBuffer calculatetoAlignmentDistance(Alignment trueAlignment) {
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < getSpectrumCount(); i++) {
			Spectrum spe = getSpectrum(i);
			double dist = 0;
			for (int j = 0; j < trueAlignment.getSequenceCount(); j++) {
				Sequence seq = trueAlignment.getSequence(j);
				for (int k = 0; k < getSpectrumLength(); k++) {
					int state = seq.getState(k);
					double d = spe.getFrequenciesAt(k)[state]; //Match
					dist += d;
				}
				dist /= getSpectrumLength();
				sb.append(dist).append("\t");
			}
			sb.append("\n");
		}
		sb.append("\n");
		return sb;
	}

	
	public static SpectrumAlignmentModel importPartialSpectrumFile(AlignmentMapping aMap,
			String partialSpectrumName) throws IOException {

		BufferedReader inFile  = new BufferedReader(new FileReader(partialSpectrumName));
		String inline;
		int length = 0;
		SpectrumAlignmentModel spectrumModel = null;
		
		while((inline = inFile.readLine())!=null){
			
			if(inline.startsWith(">")){
				String taxonName = inline.substring(1, inline.length()).trim();
				Taxon taxon = new Taxon(taxonName);
//				System.out.println(name);
				
				if(length != 0){
					double[][] freqs = new double[4][length];
					for (int i = 0; i < 4; i++) {
						inline = inFile.readLine();
//						StringTokenizer st = new StringTokenizer(inline);
						String[] result = inline.split("\\s");
						for (int j = 0; j < length; j++) {
							freqs[i][j] = Double.parseDouble(result[j]);
						}
					}
					Spectrum spectrum = new Spectrum(freqs);
					spectrum.setTaxon(taxon);
					
					int taxonIndex = spectrumModel.getTaxonIndex(taxonName);
					if(taxonIndex != -1){
						spectrumModel.removeSpectrum(taxonIndex);
						System.err.println("remove "+taxonName +"\t"+ taxonIndex);
					}
					spectrumModel.addSpectrum(spectrum);
					
					
				}
				else{
					inline = inFile.readLine();
					String[] result = inline.split("\\s");
					
					length = result.length;
					spectrumModel = new SpectrumAlignmentModel(aMap);
//					spectrumModel.aMap = aMap;
					
					
					double[][] freqs = new double[4][length];
					for (int j = 0; j < length; j++) {
						freqs[0][j] = Double.parseDouble(result[j]);
					}
				
					for (int i = 1; i < 4; i++) {
						inline = inFile.readLine();
						result = inline.split("\\s");
						for (int j = 0; j < length; j++) {
							freqs[i][j] = Double.parseDouble(result[j]);
						}
					}
					Spectrum spectrum = new Spectrum(freqs);
					spectrum.setTaxon(taxon);
					spectrumModel.addSpectrum(spectrum);
					
				}
				 
				
			}
			
				
				
		}
//		int length = 0;
//		SpectrumAlignmentModel model = new SpectrumAlignmentModel(length);
		inFile.close();
		return spectrumModel;
	}
	
}