package srp.evolution.spectrum;


import java.util.ArrayList;

import srp.evolution.OperationType;
import srp.evolution.OperationRecord;
import srp.evolution.spectrum.CategorySpectraParameter.CategoryType;
import srp.evolution.spectrum.SpectraParameter.SpectraType;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxon;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.util.NumberFormatter;



public class CategorySpectrumAlignmentModel extends AbstractSpectrumAlignmentModel  {

	
	private static final boolean DEBUG = false;
	
	private static final long serialVersionUID = 3458306765918357829L;
	private static final String MODEL_NAME = "CategorySpectrumModel";

	protected ArrayList<CategorySpectrum> categorySpectrumList;
	
	public CategorySpectrumAlignmentModel(int spectrumLength, int spectrumCount){
		this(spectrumLength, spectrumCount, CategoryType.SINGLE);
	}
	
	/**
	 * @param spectrumLength
	 * @param spectrumCount
	 * @param type
	 */
	public CategorySpectrumAlignmentModel(int spectrumLength, int spectrumCount, CategoryType type){ 
		this(spectrumLength);		
		for (int i = 0; i < spectrumCount; i++) {
			Taxon t = new Taxon(TAXON_PREFIX+i); 
			CategorySpectrum spectrum = new CategorySpectrum(spectrumLength, type);
			spectrum.setTaxon(t);
			addSpectrum(spectrum);
//			randomHaplotype(i);
		}
	}

	
	
	public CategorySpectrumAlignmentModel(int spectrumLength) {
		super(MODEL_NAME);
		this.spectrumLength = spectrumLength;
		categorySpectrumList = new ArrayList<CategorySpectrum>();
		OperationRecord = new OperationRecord();
	}

	public static CategorySpectrumAlignmentModel duplicateSpectrumAlignmentModel(CategorySpectrumAlignmentModel oldModel){
		
		CategorySpectrumAlignmentModel newSpectrumModel = new CategorySpectrumAlignmentModel(oldModel.getSpectrumLength());
		for (int i = 0; i < oldModel.getSpectrumCount(); i++) {
			CategorySpectrum spectrum = CategorySpectrum.duplicateSpectrum(oldModel.getSpectrum(i));
			newSpectrumModel.addSpectrum(spectrum);
		}
		return newSpectrumModel;
	}

	public void addSpectrum(CategorySpectrum spectrum) {
	
		for (int i = 0; i < spectrumLength; i++) {
			spectrum.setStoreSiteIndex(i);
			spectrum.storeState();
		}
		categorySpectrumList.add(spectrum);
	
	}

	public int getSpectrumCount() {
		return categorySpectrumList.size();
	}

	@Override
	public CategorySpectrum getSpectrum(int i){
		return categorySpectrumList.get(i);
	
	}

	@Override
	public String getSpectrumString(int i) {
		return categorySpectrumList.get(i).toString();
	}

	@Override
	public void removeSpectrum(int i) {
		categorySpectrumList.remove(i);
	}
	//////////////////////////////////////////////////////////////////////////
	

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
				System.out.println("StoreState in CategorySpectrumAlignment:\t"+operation);
			}
			break;
		case FULL:
			if(DEBUG){
				System.out.println("StoreState in CategorySpectrumAlignment:\t"+operation);
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
			throw new IllegalArgumentException("Unknown operation type: "+operation +"\tin"+CategorySpectrumAlignmentModel.class.getSimpleName() );
			
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
//		System.err.println("zzzzRestore CategorySpectrumAlignment: "+operation);
		switch (operation) {
		
		case NONE:
			if(DEBUG){
				System.out.println("RestoreState in CategorySpectrumAlignment:\t"+operation);
			}
			break;
		case FULL:
			if(DEBUG){
				System.out.println("RestoreState in CategorySpectrumAlignment:\t"+operation);
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