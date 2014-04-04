package srp.evolution.spectrum;


import java.util.ArrayList;

import srp.evolution.OperationType;
import srp.evolution.OperationRecord;
import srp.evolution.spectrum.SpectraParameter.SpectraType;
import dr.evolution.alignment.Alignment;
import dr.evolution.util.Taxon;
import dr.util.NumberFormatter;



public class SpectrumAlignmentModel extends AbstractSpectrumAlignmentModel  {

	
	private static final boolean DEBUG = false;
	
	private static final long serialVersionUID = 3458306765918357829L;
	private static final String MODEL_NAME = "SpectrumModel";

	protected ArrayList<Spectrum> spectrumList;
//	protected ArrayList<Spectrum> storedSpectrumList;

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
		OperationRecord = new OperationRecord();
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

	@Override
	public int getSpectrumCount() {
		return spectrumList.size();
	}

	@Override
	public Spectrum getSpectrum(int i){
		return spectrumList.get(i);
	}

	@Override
	public String getSpectrumString(int i) {
		return spectrumList.get(i).toString();
		
	}

	
	@Override
	public void removeSpectrum(int i) {
		spectrumList.remove(i);
		
	}

	

	
	@Override
	protected void acceptState() {
		//Do nothing
	}
	@Override
	protected void storeState() {

		
		OperationType operation = OperationRecord.getOperation();
		int spectrumIndex;
		int siteIndex;
		int[] siteIndexs;
		AbstractSpectrum spectrum;
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
		case COLUMN:
		case SWAP_SUBCOLUMN:
//			spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			siteIndex = OperationRecord.getSingleIndex();
			for (int i = 0; i < getSpectrumCount(); i++) {
				spectrum = getSpectrum(i);
				spectrum.setStoreSiteIndex(siteIndex);
				spectrum.storeState();
			}
			break;
		case SINGLE:
		
			spectrumIndex = OperationRecord.getSpectrumIndex();
			siteIndex = OperationRecord.getSingleIndex();
			spectrum = getSpectrum(spectrumIndex);
				spectrum.setStoreSiteIndex(siteIndex);
				spectrum.storeState();
			break;
		
		case MULTI:
		
			spectrumIndex = OperationRecord.getSpectrumIndex();
			siteIndexs = OperationRecord.getAllSiteIndexs();
			spectrum = getSpectrum(spectrumIndex);
			for (int i = 0; i < siteIndexs.length; i++) {
				spectrum.setStoreSiteIndex(siteIndexs[i]);
				spectrum.storeState();
			}
			break;
		case RECOMBINATION:

//			System.err.println("store alignment recombination");
			int[] twoSpectrums = OperationRecord.getRecombinationSpectrumIndex();
			int[] twoPositions = OperationRecord.getRecombinationPositionIndex();
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
		
		OperationType operation = OperationRecord.getOperation();
		int spectrumIndex;
		int siteIndex;
		AbstractSpectrum spectrum;
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
		case COLUMN:
		case SWAP_SUBCOLUMN:
//				spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			siteIndex = OperationRecord.getSingleIndex();
//				System.err.println(spectrumIndex +"\t"+ siteIndex +"\t"+ Arrays.toString(getSpecturmFrequencies(spectrumIndex,
//						siteIndex)));
			for (int i = 0; i < getSpectrumCount(); i++) {
				spectrum = getSpectrum(i);
				spectrum.setStoreSiteIndex(siteIndex);
				spectrum.restoreState();
			}
			break;
		case MULTI:
		
			spectrumIndex = OperationRecord.getSpectrumIndex();
			siteIndexs = OperationRecord.getAllSiteIndexs();
			spectrum = getSpectrum(spectrumIndex);
			for (int i = 0; i < siteIndexs.length; i++) {
				spectrum.setStoreSiteIndex(siteIndexs[i]);
				spectrum.restoreState();
			}
			break;
		case SINGLE:
		
			spectrumIndex = OperationRecord.getSpectrumIndex();
			siteIndex = OperationRecord.getSingleIndex();
			spectrum = getSpectrum(spectrumIndex);

				spectrum.setStoreSiteIndex(siteIndex);
				spectrum.restoreState();

			break;
			
		case RECOMBINATION:
//			System.err.println("Restore alignment recombination");
			int[] twoSpectrums = OperationRecord.getRecombinationSpectrumIndex();
			int[] twoPositions = OperationRecord.getRecombinationPositionIndex();
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



		
}