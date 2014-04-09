package srp.evolution.spectrum;

import srp.evolution.AbstractAlignmentModel;
import srp.evolution.OperationType;
import srp.evolution.OperationRecord;
import dr.evolution.alignment.SiteList;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;

//public abstract class AbstractSpectrumAlignmentModel extends AbstractModel implements TaxonList, SiteList{
public abstract class AbstractSpectrumAlignmentModel extends AbstractAlignmentModel {
//maybe only can implements TaxonList
	

	private static final long serialVersionUID = -214692337132875593L;
	protected int spectrumLength;
	public AbstractSpectrumAlignmentModel(String name) {
		super(name);
	}

	public int getSpectrumLength() {
		
		return spectrumLength;
	}
	
	public double[] getSpecturmFrequencies(int spectrumIndex, int i) {
			return getSpectrum(spectrumIndex).getFrequenciesAt(i);
	}
	
	public abstract int getSpectrumCount();
	public abstract AbstractSpectrum getSpectrum(int i);
	public abstract void removeSpectrum(int i);
	public abstract String getSpectrumString(int i);

	
	


    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

	@Override
	public int getTaxonCount() {
		return getSpectrumCount();
	}
	@Override
	public Taxon getTaxon(int taxonIndex) {
		return getSpectrum(taxonIndex).getTaxon();
	}
	




}
