package srp.evolution.spectrum;

import java.text.FieldPosition;
import java.text.NumberFormat;
import java.util.Arrays;

import srp.dr.evolution.datatype.ShortReads;
import srp.evolution.OperationRecord;
import dr.evolution.alignment.Alignment;

public class SpectrumAlignmentUtils {

	public static final int DIM = SpectraParameter.DIMENSION;
	public static final ShortReads DATA_TYPE = ShortReads.INSTANCE;

	private static final NumberFormat formatter = NumberFormat.getNumberInstance();
	
	public static double[][] compareSpectrumToTrueAlignment(AbstractSpectrumAlignmentModel spectrumModel,
			Alignment alignment, Dist dist){
		
//		Dist dist = Dist.abs;
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
			AbstractSpectrum spectrum = spectrumModel.getSpectrum(i);
			for (int j = 0; j < hapCount; j++) {
				String haplotype = alignment.getAlignedSequenceString(j);
				totalDelta[i][j] = computeDist(spectrum, haplotype, dist);
			}
			
		}
		return totalDelta;
	}
	

	public static double[][] compareSpectrum(
			AbstractSpectrumAlignmentModel spectrumModel, Dist dist) {

		int specCount = spectrumModel.getSpectrumCount();
		double[][] totalDelta = new double[specCount][specCount];
		for (int i = 0; i < specCount; i++) {
			AbstractSpectrum spectrum1 = spectrumModel.getSpectrum(i);
			for (int j = 0; j < specCount; j++) {
				AbstractSpectrum spectrum2 = spectrumModel.getSpectrum(j);
				totalDelta[i][j] = computeDist(spectrum1, spectrum2, dist);
			}
			
		}
		return totalDelta;
	}

	private static double computeDist(AbstractSpectrum spectrum1, AbstractSpectrum spectrum2, Dist dist) {
		double delta = -1;
		switch (dist) {
		case euclidean:
			delta = computeEuclideanDist(spectrum1, spectrum2);
			break;
		case abs:
			delta = computeAbsDist(spectrum1, spectrum2);
			break;
		case major:
			delta = computeMajorDist(spectrum1, spectrum2);
			break;
		default:
			break;
		}
		return delta;
	}



	private static double computeAbsDist(AbstractSpectrum spectrum1, AbstractSpectrum spectrum2) {
			int specLength = spectrum1.getLength();
			double allDist = 0;
			for (int s = 0; s < specLength; s++) {
				AbstractSpectra spectra1 = spectrum1.getSpectra(s);
				AbstractSpectra spectra2 = spectrum2.getSpectra(s);
	//			int state = getStateAtK(haplotype, s);
	//			double dist = 0;
				for (int f = 0; f < DIM; f++) {
					double delta = spectra1.getFrequency(f) - spectra2.getFrequency(f);
					allDist += Math.abs(delta);
				}
	//			allDist += Math.sqrt(dist);
			}
			allDist /= specLength;
			return allDist;
		
		}


	private static double computeEuclideanDist(AbstractSpectrum spectrum1, AbstractSpectrum spectrum2) {
			int specLength = spectrum1.getLength();
			double allDist = 0;
			for (int s = 0; s < specLength; s++) {
				AbstractSpectra spectra1 = spectrum1.getSpectra(s);
				AbstractSpectra spectra2 = spectrum2.getSpectra(s);
	//			int state = getStateAtK(haplotype, s);
				double dist = 0;
				for (int f = 0; f < DIM; f++) {
					double delta = spectra1.getFrequency(f) - spectra2.getFrequency(f);
					dist += delta*delta;
				}
				allDist += Math.sqrt(dist);
			}
			allDist /= specLength;
			return allDist;
		}


	private static double computeMajorDist(AbstractSpectrum spectrum1, AbstractSpectrum spectrum2) {
		int specLength = spectrum1.getLength();
		double allDist = 0;
		for (int s = 0; s < specLength; s++) {
			AbstractSpectra spectra1 = spectrum1.getSpectra(s);
			AbstractSpectra spectra2 = spectrum2.getSpectra(s);
//			int state = getStateAtK(haplotype, s);
//			double dist = 0;
			double max1 = 0;
			double max2 = 0;
			int maxState1 = -1;
			int maxState2 = -1;
			for (int f = 0; f < DIM; f++) {
				double delta = spectra1.getFrequency(f);
				if(delta > max1){
					max1 = delta;
					maxState1 = f;
				}
				delta = spectra2.getFrequency(f);
				if(delta > max2){
					max2 = delta;
					maxState2 = f;
				}
			}
			
			if(maxState1 != maxState2){
				allDist++;
			}
//			allDist += Math.sqrt(dist);
		}
		allDist /= specLength;
		return allDist;
	
	}

	private static double computeDist(AbstractSpectrum spectrum, String haplotype, Dist dist) {
		double delta = -1;
		switch (dist) {
		case euclidean:
			delta = computeEuclideanDist(spectrum, haplotype);
			break;
		case abs:
			delta = computeAbsDist(spectrum, haplotype);
			break;
		case major:
			delta = computeMajorDist(spectrum, haplotype);
			break;
		default:
			break;
		}
		return delta;
	}
	
	private static double computeAbsDist(AbstractSpectrum spectrum, String haplotype) {
		int specLength = spectrum.getLength();
		double dist = 0;
		for (int s = 0; s < specLength; s++) {
			AbstractSpectra spectra = spectrum.getSpectra(s);
			int state = getStateAtK(haplotype, s);
			for (int f = 0; f < DIM; f++) {
				double delta = spectra.getFrequency(f);
				if(f==state){
					delta = 1-delta;
				}
//				dist += Math.abs(delta);
				dist += delta;
			}
		}
		dist /= specLength;
		return dist;
	}
	
	private static double computeEuclideanDist(AbstractSpectrum spectrum, String haplotype) {
		int specLength = spectrum.getLength();
		double allDist = 0;
		for (int s = 0; s < specLength; s++) {
			AbstractSpectra spectra = spectrum.getSpectra(s);
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
//			System.out.println("D: "+dist +"\t"+allDist +"\n");
		}
//		System.exit(-1);
		allDist /= specLength;
		return allDist;
	}


	private static double computeMajorDist(AbstractSpectrum spectrum, String haplotype) {
			int specLength = spectrum.getLength();
			double dist = 0;
			for (int s = 0; s < specLength; s++) {
				AbstractSpectra spectra = spectrum.getSpectra(s);
				int state = getStateAtK(haplotype, s);
				double max = 0;
				int maxState = -1;
				for (int f = 0; f < DIM; f++) {
					double delta = spectra.getFrequency(f);
					if(delta > max){
						max = delta;
						maxState = f;
					}
	//				if(f==state){
	//					delta = 1-delta;
	//				}
	//				dist += Math.abs(delta);
	//				dist += delta;
				}
				if(state!=maxState){
					dist++;
				}
			}
			dist /= specLength;
			return dist;
		}


	private static int getStateAtK(String fullSrp, int k) {
		char srpChar = fullSrp.charAt(k);
		int state = DATA_TYPE.getState(srpChar);
	
		return state;
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
		sb.append("\n");
		return sb.toString();
	}


	public static void compareTwoSpectrumModel(
			AbstractSpectrumAlignmentModel spectrumModel,
			AbstractSpectrumAlignmentModel spectrumModel2) {
	
	
		for (int i = 0; i < spectrumModel.getSpectrumCount(); i++) {
			AbstractSpectrum spectrum = spectrumModel.getSpectrum(i);
			AbstractSpectrum newSpectrum = spectrumModel2.getSpectrum(i);
			for (int j = 0; j < spectrum.getLength(); j++) {
				double[] frequencies = spectrum.getFrequenciesAt(j);
				double[] newFrequencies = newSpectrum.getFrequenciesAt(j);
				for (int k = 0; k < frequencies.length; k++) {
					if(frequencies[k] != newFrequencies[k]){
						System.err.println("DIFFMODEL\t"+i +"\t"+ j +"\t"+ k +"\t"+ Arrays.toString(frequencies) +"\t"+ Arrays.toString(newSpectrum.getFrequenciesAt(j)));
					}
				}
			}
			
		}
		OperationRecord record = spectrumModel.getOperationRecord();
		int spec = record.getSpectrumIndex();
		int site = record.getSingleIndex();
		System.out.println("compare two models");
		System.out.println(Arrays.toString(spectrumModel.getSpectrum(spec).getFrequenciesAt(
				site)));
		System.out.println(Arrays.toString(spectrumModel2.getSpectrum(spec).getFrequenciesAt(
				site)));
		
	}
	
	public enum Dist{
		abs,
		euclidean,
		major,
	}

	public static String CreateMajorAlignment(
			AbstractSpectrumAlignmentModel spectrumModel) {

		int specLength = spectrumModel.getSpectrumLength();
		int specCount = spectrumModel.getSpectrumCount();
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < specCount; i++) {

			AbstractSpectrum spectrum = spectrumModel.getSpectrum(i);
			for (int s = 0; s < specLength; s++) {
				AbstractSpectra spectra = spectrum.getSpectra(s);
				double max = 0;
				int maxState = -1;

				for (int f = 0; f < DIM; f++) {
					double delta = spectra.getFrequency(f);
					if(delta > max){
						max = delta;
						maxState = f;
					}
				}
				sb.append( ShortReads.NUCLEOTIDE_CHARS[maxState] );
			}
			sb.append("\n");
		}
		sb.append("\n");
		
		return sb.toString();
	}
}
