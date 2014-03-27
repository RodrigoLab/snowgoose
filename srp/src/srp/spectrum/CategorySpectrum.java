package srp.spectrum;

import java.util.ArrayList;

import srp.spectrum.CategorySpectraParameter.CategoryType;

//public class Spectrum implements Identifiable, Attributable{// extends AbstractModel{
public class CategorySpectrum extends AbstractSpectrum {

	private static final long serialVersionUID = -728370884996776301L;

	private ArrayList<CategorySpectraParameter> spectrum;

	private CategorySpectraParameter[] spectrumArray;


	public CategorySpectrum(int length) {
		this(length, CategoryType.SINGLE);
	}

	public CategorySpectrum(int length, CategoryType type) {

		super("CategorySpectrum");
		this.spectrumLength = length;
		spectrumArray = new CategorySpectraParameter[length];
		spectrum = new ArrayList<CategorySpectraParameter>();
		for (int i = 0; i < length; i++) {
			CategorySpectraParameter spectra = new CategorySpectraParameter(
					type);
			addVariable(spectra);
			spectrum.add(spectra);
			spectrumArray[i] = spectra;
		}

	}

	public void resetCategorySpectra(int site, int cat) {
		getSpectra(site).setCategory(cat);
	}
	
	@Override
	public CategorySpectraParameter getSpectra(int site) {
		return spectrum.get(site);
	}

	public CategorySpectraParameter getSpectraArray(int site) {
		return spectrumArray[site];
	}

}
