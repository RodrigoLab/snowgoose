package srp.spectrum;


import java.util.ArrayList;

import srp.dr.evolution.datatype.ShortReads;

import srp.spectrum.SpectraParameter.SpectraType;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.Nucleotides;
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

//	private static final DataType DATA_TYPE = Nucleotides.INSTANCE;
	public static final char GAP = '-';
	public static final String TAXON_PREFIX = "taxa_";

//	int haplotypesCount;


	private boolean isEdit;//TODO utilise this!!
	
	private SpectrumOperationRecord spectrumOperationRecord;
	

	public SpectrumAlignmentModel(int spectrumLength, int spectrumCount){
		this(spectrumLength, spectrumCount, SpectraType.RANDOM);
	}
	
	/**
	 * type 0=Equal. 
	 * 		1=[1 0 0 0].
	 * 		2=[Random]. 
	 * @param spectrumLength
	 * @param spectrumCount
	 * @param type
	 */
	public SpectrumAlignmentModel(int spectrumLength, int spectrumCount, SpectraType type) {
		this(spectrumLength);		
		for (int i = 0; i < spectrumCount; i++) {
			Taxon t = new Taxon(TAXON_PREFIX+i); 
			Spectrum spectrum = new Spectrum(spectrumLength, type);
			spectrum.setTaxon(t);
			addSpectrum(spectrum);
//			randomHaplotype(i);
		}
	}

	
	
	public SpectrumAlignmentModel(Alignment trueAlignment) {
		this(trueAlignment.getSiteCount());

		for (int i = 0; i < trueAlignment.getSequenceCount(); i++) {
			Spectrum spectrum = new Spectrum(trueAlignment.getSequence(i));
			Taxon t = new Taxon(TAXON_PREFIX+i); 
//			Spectrum spectrum = new Spectrum(spectrumLength);
			spectrum.setTaxon(t);
			addSpectrum(spectrum);
		}
	}
	
	
	
	public SpectrumAlignmentModel(int spectrumLength) {
		super(MODEL_NAME);
		this.spectrumLength = spectrumLength;
		
		spectrumList = new ArrayList<Spectrum>();
//		storedSpectrumList = new ArrayList<Spectrum>();
//		setDataType(ShortReads.INSTANCE);
		spectrumOperationRecord = new SpectrumOperationRecord();
	}

	public static SpectrumAlignmentModel duplicateSpectrumAlignmentModel(SpectrumAlignmentModel oldModel){
		
		SpectrumAlignmentModel newSpectrumModel = new SpectrumAlignmentModel(oldModel.getSpectrumLength());

		for (int i = 0; i < oldModel.getSpectrumCount(); i++) {
			Spectrum spectrum = Spectrum.duplicateSpectrum(oldModel.getSpectrum(i));
			newSpectrumModel.addSpectrum(spectrum);
		}
		
		
		return newSpectrumModel;
	}
	
	public void addSpectrum(Spectrum spectrum) {

		for (int i = 0; i < spectrumLength; i++) {
			spectrum.setStoreSiteIndex(i);
			spectrum.storeState();
		}
		spectrumList.add(spectrum);

	}

	public void resetSpectrumOperation(){
		spectrumOperationRecord.setOperation(SpectrumOperation.FULL);
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


	//////////////////////////////////////////////////////////////////////////
	
	
	
	public void removeSpectrum(int i) {
		spectrumList.remove(i);
		
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
	public void setSpectrumOperationRecord(SpectrumOperation op, int spectrumIndex,
			int siteIndex) {
		spectrumOperationRecord.setRecord(op, spectrumIndex, siteIndex);
	}
	public void setSpectrumOperationRecord(SpectrumOperation op, int spectrumIndex,
			int siteIndex, double delta) {
		spectrumOperationRecord.setRecord(op, spectrumIndex, siteIndex);
	}
	public void setSpectrumOperationRecord(SpectrumOperation op,
			int spectrumIndex, int[] siteIndexs, double... delta) {
		spectrumOperationRecord.setRecord(op, spectrumIndex, siteIndexs, delta);
	
	}

	//	public void setSpectrumOperationRecord(SpectrumOperation op,
	//			int spectrumIndex, int[] siteIndexs) {
	//		//recombination
	//		spectrumOperationRecord.setRecord(op, spectrumIndex, siteIndexs);
	//
	//	}
	
		public void setSpectrumOperationRecord(SpectrumOperation op,
				int[] spectrumIndexs, int siteIndex) {
			//subcolumn
			spectrumOperationRecord.setRecord(op, spectrumIndexs, siteIndex);
			
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
	@SuppressWarnings("rawtypes")
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
				System.out.println("StoreState in SpectrumAlignment:\t"+operation);
			}
			break;
		case FULL:
			if(DEBUG){
				System.out.println("StoreState in SpectrumAlignment:\t"+operation);
			}
			for (int i = 0; i < getSpectrumCount(); i++) {
				spectrum = getSpectrum(i);
				for (int s = 0; s < getSpectrumLength(); s++) {
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
			siteIndex = spectrumOperationRecord.getSingleIndex();
			for (int i = 0; i < getSpectrumCount(); i++) {
				spectrum = getSpectrum(i);
				spectrum.setStoreSiteIndex(siteIndex);
				spectrum.storeState();
			}
			break;
		case DELTA_SINGLE:
		case SWAP_SINGLE:
			spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			siteIndex = spectrumOperationRecord.getSingleIndex();
			spectrum = getSpectrum(spectrumIndex);
				spectrum.setStoreSiteIndex(siteIndex);
				spectrum.storeState();
			break;
		
		case DELTA_MULTI:
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
				for (int s = 0; s < getSpectrumLength(); s++) {
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
			siteIndex = spectrumOperationRecord.getSingleIndex();
//				System.err.println(spectrumIndex +"\t"+ siteIndex +"\t"+ Arrays.toString(getSpecturmFrequencies(spectrumIndex,
//						siteIndex)));
			for (int i = 0; i < getSpectrumCount(); i++) {
				spectrum = getSpectrum(i);
				spectrum.setStoreSiteIndex(siteIndex);
				spectrum.restoreState();
			}
			break;
		case DELTA_MULTI:
		case SWAP_MULTI:
			spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			siteIndexs = spectrumOperationRecord.getAllSiteIndexs();
			spectrum = getSpectrum(spectrumIndex);
			for (int i = 0; i < siteIndexs.length; i++) {
				spectrum.setStoreSiteIndex(siteIndexs[i]);
				spectrum.restoreState();
			}
			break;
		case DELTA_SINGLE:
		case SWAP_SINGLE:
			spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			siteIndex = spectrumOperationRecord.getSingleIndex();
			spectrum = getSpectrum(spectrumIndex);

				spectrum.setStoreSiteIndex(siteIndex);
				spectrum.restoreState();

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

    

	//REMOVE after
	private static final int NUCLEOTIDE_STATES[] = Nucleotides.NUCLEOTIDE_STATES;
	private static final char[] VALID_CHARS = initValidChars4();
	private static final int INDEX_OF_LAST_VALID_CHARS = VALID_CHARS.length-1;
	private static char[] initValidChars4() {
		char[] validChar = new char[4];
		System.arraycopy(Nucleotides.NUCLEOTIDE_CHARS, 0, validChar, 0, validChar.length);
		return validChar;
	}
		
}