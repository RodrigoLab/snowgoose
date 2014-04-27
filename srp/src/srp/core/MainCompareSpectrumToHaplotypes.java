package srp.core;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;

import srp.evolution.haplotypes.SPSDist;
import srp.evolution.spectrum.AbstractSpectrumAlignmentModel;
import srp.evolution.spectrum.SpectrumAlignmentUtils;
import srp.evolution.spectrum.SpectrumAlignmentUtils.Dist;
import dr.evolution.alignment.Alignment;

public class MainCompareSpectrumToHaplotypes {

	public static void main(String[] args) throws Exception {

		int runIndex;
		int noOfTrueHaplotype;
		int noOfRecoveredHaplotype;

		String dataOnlyDir = "/home/sw167/workspaceSrp/Simulations/dataOnly/H7/";
//		String resultDir = "/home/sw167/workspaceSrp/Simulations/srpResults/Run0314_betaMean/";
//		String resultDir = "/home/sw167/workspaceSrp/Simulations/srpResults/Run0314_flat/";
		String resultDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/testData/";
//		String resultDir = "/home/sw167/workspaceSrp/snowgoose/jar/";
		
		noOfTrueHaplotype = 7;
		noOfRecoveredHaplotype = 7;

		int i = 54;
//		for (int i = 51; i < 60; i++)
//		for (int i = 70; i < 80; i++)
		{
			runIndex = i;
			String hapRunIndex = "H" + noOfTrueHaplotype + "_" + runIndex;
			
			String dataDir = resultDir + hapRunIndex+ "/";
			String partialSpectrumName = dataDir + "FullTree_" + hapRunIndex + ".haplotype";
			String compareResultFile =  dataDir + "FullTree_" + hapRunIndex + ".dist";
			
			String trueHaplotypeFile = dataOnlyDir + hapRunIndex + "/" + hapRunIndex + "_Srp_fullHaplotype.fasta";

			System.out.println(partialSpectrumName);
			System.out.println(trueHaplotypeFile);
			// String partialSpectrumName = hapRunIndex+".haplatypepartial";
			// String partialTreeName = "FullTree_"+hapRunIndex+".treespartial";

			DataImporter dataImporter = new DataImporter("");

			Alignment trueAlignment = dataImporter
					.importAlignment(trueHaplotypeFile);

			AbstractSpectrumAlignmentModel spectrumModel = dataImporter
					.importPartialSpectrumFile(partialSpectrumName);
			
//			 loggers[3] = new SpectrumLogger(spectrumModel, trueAlignment,
//			 logHaplotypeName, logInterval*totalSamples/10);
			PrintWriter out = new PrintWriter(new BufferedWriter(
					new FileWriter(compareResultFile)));
			
			double[][] dist = SpectrumAlignmentUtils
					.compareSpectrumToTrueAlignment(spectrumModel,
							trueAlignment, Dist.abs);
			String s = SpectrumAlignmentUtils.formatter(dist);
			out.println("Abs");
			out.println(s);


			dist = SpectrumAlignmentUtils.compareSpectrumToTrueAlignment(
					spectrumModel, trueAlignment, Dist.euclidean);
			s = SpectrumAlignmentUtils.formatter(dist);
			out.println("Euclidean");
			out.println(s);

			dist = SpectrumAlignmentUtils.compareSpectrumToTrueAlignment(
					spectrumModel, trueAlignment, Dist.major);
			s = SpectrumAlignmentUtils.formatter(dist);
			out.println("Major");
			out.println(s);

			dist = SpectrumAlignmentUtils.compareSpectrum(spectrumModel,
					Dist.euclidean);
			s = SpectrumAlignmentUtils.formatter(dist);
			out.println("Within Euclidean");
			out.println(s);

			dist = SpectrumAlignmentUtils.compareSpectrum(spectrumModel,
					Dist.major);
			s = SpectrumAlignmentUtils.formatter(dist);
			out.println("Within major");
			out.println(s);

			String alignmentString = SpectrumAlignmentUtils.CreateMajorAlignment(spectrumModel);
			out.println(alignmentString);
			
			s = SPSDist.calculeteSPSArrayString(trueAlignment,
					trueAlignment);
			out.println("Within haplotype");
			out.println(s);

			out.println("================================\n");
			out.close();
		}

	}
}
