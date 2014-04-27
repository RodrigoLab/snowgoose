package srp.evolution.spectrum;

import java.util.ArrayList;

import srp.evolution.spectrum.SpectraParameter.SpectraType;
import dr.evolution.sequence.Sequence;

//public class Spectrum implements Identifiable, Attributable{// extends AbstractModel{
public class Spectrum extends AbstractSpectrum {

	private static final long serialVersionUID = -728370884996776301L;

	private ArrayList<SpectraParameter> spectrum;
	private SpectraParameter[] spectrumArray;


	public Spectrum(int length) {
		this(length, SpectraType.EQUAL);
	}

	/**
	 * type EQUAL = Equal. ZERO_ONE= [1 0 0 0]. RANDOM = [Random].
	 * 
	 * @param length
	 * @param type
	 */
	public Spectrum(int length, SpectraType type) {

		super("Spectrum");
		this.spectrumLength = length;
		spectrumArray = new SpectraParameter[length];
		spectrum = new ArrayList<SpectraParameter>();
		for (int i = 0; i < length; i++) {
			SpectraParameter spectra = new SpectraParameter(type);
			addVariable(spectra);
			spectrum.add(spectra);
			spectrumArray[i] = spectra;
		}

	}

	public Spectrum(double[][] freqs) {
		this(freqs[0].length);
		for (int s = 0; s < spectrum.size(); s++) {
			SpectraParameter spectra = spectrum.get(s);
			for (int f = 0; f < 4; f++) {
				spectra.setFrequency(f, freqs[f][s]);
			}

		}
	}

	public Spectrum(Sequence sequence) {
		this(sequence.getLength());

		for (int i = 0; i < sequence.getLength(); i++) {
			SpectraParameter spectra = getSpectra(i);
			char c = sequence.getChar(i);
			int state = DATA_TYPE.getState(c);
			for (int j = 0; j < STATE_COUNT; j++) {
				if (j == state) {
					spectra.setFrequency(j, 1);
				} else {
					spectra.setFrequency(j, 0);
				}
			}
		}

	}

	

	@Override
	public SpectraParameter getSpectra(int site) {
		return spectrum.get(site);
	}

	
	public void resetSpectra(int site, double[] values) {
		SpectraParameter spectra = getSpectra(site);//.setFrequenciesQuietly(values);
		for (int i = 0; i < SpectraParameter.DIMENSION; i++) {
			spectra.setFrequency(i, values[i]);
		}
	}

	public SpectraParameter getSpectraArray(int site) {
		return spectrumArray[site];
	}

}
