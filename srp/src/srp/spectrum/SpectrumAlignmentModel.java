package srp.spectrum;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import srp.haplotypes.AlignmentMapping;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.util.NumberFormatter;



public class SpectrumAlignmentModel extends AbstractSpectrumAlignmentModel  {

	
	private static final boolean DEBUG = false;
	
	private static final long serialVersionUID = 3458306765918357829L;
	private static final String MODEL_NAME = "SpectrumModel";

	
	private static final int NUCLEOTIDE_STATES[] = Nucleotides.NUCLEOTIDE_STATES;
	private static final char[] VALID_CHARS = initValidChars4();
	private static final int INDEX_OF_LAST_VALID_CHARS = VALID_CHARS.length-1;
	
//	private static final DataType DATA_TYPE = Nucleotides.INSTANCE;
	public static final char GAP = '-';
	public static final String TAXON_PREFIX = "taxa_";

//	int haplotypesCount;
	private AlignmentMapping aMap;

	private boolean isEdit;
	
	private SpectrumOperationRecord spectrumOperationRecord;
	
//	public long time = 0;

	

	
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

	@Deprecated
	public SpectrumAlignmentModel(Alignment shortReads, int hapCount) {
		this(new AlignmentMapping(shortReads), hapCount);
	}
	@Deprecated
	public SpectrumAlignmentModel(AlignmentMapping aMap, int hapCount) {
		this(aMap, hapCount, 1);		
	}
	/**
	 * type 0=Equal. 
	 * 		1=[1 0 0 0].
	 * 		2=[Random]. 
	 * @param aMap
	 * @param hapCount
	 * @param type
	 */
	public SpectrumAlignmentModel(AlignmentMapping aMap, int hapCount, int type) {
		this(aMap);		
		for (int i = 0; i < hapCount; i++) {
			Taxon t = new Taxon(TAXON_PREFIX+i); 
			Spectrum spectrum = new Spectrum(spectrumLength, type);
			spectrum.setTaxon(t);
			addSpectrum(spectrum);
//			randomHaplotype(i);
		}
	}

	public SpectrumAlignmentModel(AlignmentMapping aMap, int hapCount, Type type) {
		this(aMap, hapCount, type.getCode());
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


	public void resetSpectrumOperation(){
		spectrumOperationRecord.setOperation(SpectrumOperation.FULL);
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
	public void setSpectrumOperationRecord(SpectrumOperation op,
			int spectrumIndex, int[] siteIndexs) {
		//recombination
		spectrumOperationRecord.setRecord(op, spectrumIndex, siteIndexs);

	}

	public void setSpectrumOperationRecord(SpectrumOperation op,
			int[] spectrumIndexs, int siteIndex) {
		//subcolumn
		spectrumOperationRecord.setRecord(op, spectrumIndexs, siteIndex);
		
	}


	public SpectrumOperationRecord getSpectrumOperationRecord() {
		return spectrumOperationRecord;
	}

	public SpectrumOperation getSpectrumOperation() {
		return spectrumOperationRecord.getOperation();
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
			if(DEBUG){
				System.out.println("StoreState in SpectrumAlignment:\ts"+operation);
			}
			break;
		case FULL:
			if(DEBUG){
				System.out.println("StoreState in SpectrumAlignment:\t"+operation);
			}
			for (int i = 0; i < getSpectrumCount(); i++) {
				spectrum = getSpectrum(i);
				for (int s = 0; s < spectrumLength; s++) {
					spectrum.setStoreSiteIndex(s);
					spectrum.storeState();
				}
			}
			break;
//		case SINGLE_DELTA:
//			spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
//			siteIndex = spectrumOperationRecord.getAllSiteIndexs()[0];
//
//			spectrum = getSpectrum(spectrumIndex);
//			spectrum.setStoreSiteIndex(siteIndex);
//			spectrum.storeState();
//			break;
		case DELTA_COLUMN:
		case SWAP_COLUMN:
		case SWAP_SUBCOLUMN:
//			spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			siteIndex = spectrumOperationRecord.getColumnIndex();
			for (int i = 0; i < getSpectrumCount(); i++) {
				spectrum = getSpectrum(i);
				spectrum.setStoreSiteIndex(siteIndex);
				spectrum.storeState();
			}
			break;
		case DELTA_SINGLE:
		case DELTA_MULTI:
		case SWAP_SINGLE:
		case SWAP_MULTI:
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
				for (int s = twoPositions[0]; s < twoPositions[1]; s++) {
					spectrum.setStoreSiteIndex(s);
					spectrum.storeState();
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
		
		SpectrumOperation operation = spectrumOperationRecord.getOperation();
		int spectrumIndex;
		int siteIndex;
		Spectrum spectrum;
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
			for (int i = 0; i < getSpectrumCount(); i++) {
				spectrum = getSpectrum(i);
				for (int s = 0; s < spectrumLength; s++) {
					spectrum.setStoreSiteIndex(s);
					spectrum.restoreState();
				}
			}
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
		case DELTA_COLUMN:
		case SWAP_COLUMN:
		case SWAP_SUBCOLUMN:
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
		case DELTA_SINGLE:
		case DELTA_MULTI:
		case SWAP_SINGLE:
		case SWAP_MULTI:
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
				for (int s = twoPositions[0]; s < twoPositions[1]; s++) {
					spectrum.setStoreSiteIndex(s);
					spectrum.restoreState();
				}
			}
			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: " + operation);

		}
//		long time2 = System.currentTimeMillis();
//		time += (time2-time1);
	}
	public void startSpectrumOperation(){
//		System.err.println("\n!!!startSpectrumOperation");
		isEdit = true;
	}
	//TODO: Implement check system here
	public void endSpectrumOperation(){
		isEdit = false;
//		System.err.println("!!!!!!!endSpectrumOperation, fireModelCHanged");
		fireModelChanged();
	}
//////////////////////////////////////////////////////////////////////////

	

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


	public enum Type{
		
		EQUAL(0),
		ZERO_ONE(1),
		RANDOM(2);
		int type;
		private Type(int t){
			type = t;
		}
		public int getCode(){
			return type;
		}
	}
}