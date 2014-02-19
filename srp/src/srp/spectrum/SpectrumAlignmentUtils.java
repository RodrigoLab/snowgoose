package srp.spectrum;

import java.text.FieldPosition;
import java.text.NumberFormat;
import java.util.Arrays;

import srp.core.DataImporter;
import srp.dr.evolution.datatype.ShortReads;
import srp.haplotypes.AlignmentMapping;
import srp.spectrum.SpectrumAlignmentModel.SpectrumType;
import dr.evolution.alignment.Alignment;

public class SpectrumAlignmentUtils {

	public static final int DIM = SpectraParameter.DIMENSION;
	
	public static ShortReads dataType = ShortReads.INSTANCE;
	private final static NumberFormat formatter = NumberFormat.getNumberInstance();
	
	public static void main(String[] args) throws Exception {
		String dataDir;

		int runIndex = 50;
		int noOfTrueHaplotype = 7;
		String hapRunIndex = "H"+noOfTrueHaplotype+"_"+runIndex;
		dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/testData/"+hapRunIndex+"/";
		
		
		String shortReadFile = hapRunIndex +"_Srp.fasta";
		String trueHaplotypeFile = hapRunIndex +"_Srp_fullHaplotype.fasta";
//		String trueHaplotypeFile = "H7_50_Srp_fullHaplotype.fasta"

		String storeSpectrumDir = dataDir+"specturm_betaMean_boundry_good/";
		String prefix = storeSpectrumDir+"FullTree_"+hapRunIndex;
		String logSpectrumName = prefix+".haplatype";

//		System.exit(1);

		DataImporter dataImporter = new DataImporter(dataDir);
		Alignment trueAlignment = dataImporter.importAlignment(trueHaplotypeFile);
		Alignment shortReads = dataImporter.importShortReads(shortReadFile);
		AlignmentMapping aMap = new AlignmentMapping(shortReads);
		
		dataImporter = new DataImporter(storeSpectrumDir);
		SpectrumAlignmentModel spectrumModel = SpectrumAlignmentModel.importPartialSpectrumFile(aMap, logSpectrumName );
//			
		compareSpectrumToTrueAlignment(spectrumModel, trueAlignment);


	}
	
	public static String compareSpectrumToTrueAlignment(SpectrumAlignmentModel spectrumModel,
			Alignment alignment){
		
		StringBuffer sb = new StringBuffer();
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
		
		System.out.println(sb.toString());
		
		return sb.toString();
	}
	
	private static double computeEuclideanDist(Spectrum spectrum, String haplotype) {
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
				dist += delta*delta;

			}
		}
		dist = Math.sqrt(dist);
		return dist;
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
		int state = dataType.getState(srpChar);
	
		return state;
	}
}
