package srp.spectrum;

import java.text.FieldPosition;
import java.text.NumberFormat;
import java.util.Arrays;

import javax.swing.text.TabableView;

import srp.core.DataImporter;
import srp.dr.evolution.datatype.ShortReads;
import srp.shortreads.AlignmentMapping;
import srp.spectrum.SpectraParameter.SpectraType;
import dr.evolution.alignment.Alignment;

public class SpectrumAlignmentUtils {

	public static final int DIM = SpectraParameter.DIMENSION;
	public static final ShortReads DATA_TYPE = ShortReads.INSTANCE;

	private static final NumberFormat formatter = NumberFormat.getNumberInstance();
	
	public static double[][] compareSpectrumToTrueAlignment(SpectrumAlignmentModel spectrumModel,
			Alignment alignment){
		
		
		int hapCount = alignment.getSequenceCount();
		int hapLength = alignment.getSiteCount();
		int specCount = spectrumModel.getSpectrumCount();
		int specLength = spectrumModel.getSpectrumLength();

		if(specLength != alignment.getSiteCount()){
			System.err.println("Incompariable alignments lenght: spectrum:" + specLength
					+ " and alignment:" + alignment.getSiteCount());
			System.exit(-1);
		}

		if (specCount != hapCount) {
			System.err.println("TODO: Implement comparing different size later");
			throw new IllegalArgumentException("Not yet implemented");
		}
		
		double[][] totalDelta = new double[specCount][hapCount];
		for (int i = 0; i < specCount; i++) {
			Spectrum spectrum = spectrumModel.getSpectrum(i);
			for (int j = 0; j < hapCount; j++) {
				String haplotype = alignment.getAlignedSequenceString(j);
				totalDelta[i][j] = computeEuclideanDist(spectrum, haplotype);
			}
			
		}
		return totalDelta;
	}
	public static String formatter(double[][] totalDelta){
	
		StringBuffer sb = new StringBuffer();
		
		formatter.setMaximumFractionDigits(3);
		formatter.setMinimumFractionDigits(3);
		FieldPosition fp = new FieldPosition(12);
		for (int i = 0; i < totalDelta.length; i++) {
			for (int j = 0; j < totalDelta[i].length; j++) {
				formatter.format(totalDelta[i][j], sb, fp);
				sb.append("\t");
			}
			sb.append("\n");
		}
		return sb.toString();
	}
	
	private static double computeEuclideanDist(Spectrum spectrum, String haplotype) {
		int specLength = spectrum.getLength();
		double allDist = 0;
		for (int s = 0; s < specLength; s++) {
			SpectraParameter spectra = spectrum.getSpectra(s);
			int state = getStateAtK(haplotype, s);
			double dist = 0;
			for (int f = 0; f < DIM; f++) {
				double delta = spectra.getFrequency(f);
				if(f==state){
					delta -= 1;
				}
				dist += delta*delta;

			}
			allDist += Math.sqrt(dist);
		}
		
		return allDist;
	}
	
	private static double computeAbsDist(Spectrum spectrum, String haplotype) {
		int specLength = spectrum.getLength();
		double dist = 0;
		for (int s = 0; s < specLength; s++) {
			SpectraParameter spectra = spectrum.getSpectra(s);
			int state = getStateAtK(haplotype, s);
			for (int f = 0; f < DIM; f++) {
				double delta = spectra.getFrequency(f);
				if(f==state){
					delta -= 1;
				}
				dist += Math.abs(delta);
			}
		}
		return dist;
	}
	private static int getStateAtK(String fullSrp, int k) {
		char srpChar = fullSrp.charAt(k);
		int state = DATA_TYPE.getState(srpChar);
	
		return state;
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
}
