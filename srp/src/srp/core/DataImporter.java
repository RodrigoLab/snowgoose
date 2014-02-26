package srp.core;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import srp.dr.evolution.datatype.ShortReads;
import srp.dr.evolution.io.ShortReadImporter;
import srp.haplotypes.AlignmentMapping;
import srp.spectrum.SpectraParameter;
import srp.spectrum.Spectrum;
import srp.spectrum.SpectrumAlignmentModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.FastaImporter;
import dr.evolution.io.NewickImporter;
//import dr.evolution.io.Importer.ImportException;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.SequenceList;
import dr.evolution.sequence.Sequences;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;

public class DataImporter {

//	private SimpleAlignment alignment;
//	private FlexibleTree tree;
	
	private static final double MIN_FREQ = SpectraParameter.MIN_FREQ;
	private String dataDir;
	
	public DataImporter(String dataDir){
		if(dataDir.lastIndexOf(File.separator)!= dataDir.length()){
			dataDir += File.separator;
		}
		this.dataDir = dataDir;
	}
	
//	public DataImporter() {
//		this.dataDir = "";
//	}
	public Alignment importShortReads(String fileName) throws Exception{
		
		Alignment alignment = importShortReads(dataDir, fileName);
		return alignment;
	}
	
	public Alignment importAlignment(String fileName) throws Exception{
		
		Alignment alignment = importAlignment(dataDir, fileName);
		return alignment;
	}
	
	public Tree importTree(String fileName) throws Exception{
		
		Tree tree =  importTree(dataDir, fileName);
		return tree;
	}

	public Sequences importSequence(String fileName) throws Exception{
		
		Sequences seqs =  (Sequences) importSequences(dataDir, fileName);
		return seqs;
	}
	
	public Sequence importRefSeq(String fileName) throws Exception{
		Sequence seq = importRefSeq(dataDir, fileName);
		return seq;
	}
	
	public SpectrumAlignmentModel importPartialSpectrumFile(String partialSpectrumName) throws IOException {
		
		return importPartialSpectrumFile(dataDir, partialSpectrumName);
	}

	public static SpectrumAlignmentModel importPartialSpectrumFile(
			String dataDir, String partialSpectrumName) throws IOException {
	
		BufferedReader inFile  = new BufferedReader(new FileReader(dataDir+partialSpectrumName));
		String inline;
		int length = 0;
		SpectrumAlignmentModel spectrumModel = null;
		
		while((inline = inFile.readLine())!=null){
			
			if(inline.startsWith(">")){
				String taxonName = inline.substring(1, inline.length()).trim();
				Taxon taxon = new Taxon(taxonName);
//				System.out.println(name);
				
				if(length != 0){
					double[][] freqs = new double[4][length];
					for (int i = 0; i < 4; i++) {
						inline = inFile.readLine();
//						StringTokenizer st = new StringTokenizer(inline);
						String[] result = inline.split("\\s");
						for (int j = 0; j < length; j++) {
							freqs[i][j] = Double.parseDouble(result[j]);
							if(freqs[i][j]< MIN_FREQ){
								freqs[i][j]= MIN_FREQ;
							}
						}
					}
					Spectrum spectrum = new Spectrum(freqs);
					spectrum.setTaxon(taxon);
					
					int taxonIndex = spectrumModel.getTaxonIndex(taxonName);
					if(taxonIndex != -1){
						spectrumModel.removeSpectrum(taxonIndex);
						System.err.println("remove "+taxonName +"\t"+ taxonIndex);
					}
					spectrumModel.addSpectrum(spectrum);
					
					
				}
				else{
					inline = inFile.readLine();
					String[] result = inline.split("\\s");
					
					length = result.length;
					spectrumModel = new SpectrumAlignmentModel(length);
					
					double[][] freqs = new double[4][length];
					for (int j = 0; j < length; j++) {
						freqs[0][j] = Double.parseDouble(result[j]);
						if(freqs[0][j]< MIN_FREQ){
							freqs[0][j]=MIN_FREQ;
						}
					}
				
					for (int i = 1; i < 4; i++) {
						inline = inFile.readLine();
						result = inline.split("\\s");
						for (int j = 0; j < length; j++) {
							freqs[i][j] = Double.parseDouble(result[j]);
							if(freqs[i][j]< MIN_FREQ){
								freqs[i][j] = MIN_FREQ;
							}
						}
					}
					Spectrum spectrum = new Spectrum(freqs);
					spectrum.setTaxon(taxon);
					spectrumModel.addSpectrum(spectrum);
					
				}
				 
				
			}
			
				
				
		}
//		int length = 0;
//		SpectrumAlignmentModel model = new SpectrumAlignmentModel(length);
		inFile.close();
		return spectrumModel;
	}

	public static Sequence importRefSeq(String dataDir, String fileName) throws Exception{

		
		String filePath = dataDir+fileName;
		
		BufferedReader in = new BufferedReader(new FileReader(filePath));
		FastaImporter importer = new FastaImporter(in, Nucleotides.INSTANCE);
		SequenceList seqList = importer.importSequences();
		Sequence seq = seqList.getSequence(0);
		if(seqList.getSequenceCount()>1){
			System.err.println("Multiple ref seq\t"+seqList.getSequenceCount());
		}
	
		return seq;
	}

	public static SequenceList importSequences(String dataDir, String fileName) throws Exception {

		
		String filePath = dataDir+fileName;
		
		BufferedReader in = new BufferedReader(new FileReader(filePath));
		FastaImporter importer = new FastaImporter(in, Nucleotides.INSTANCE);
		SequenceList seqs = importer.importSequences();

		return seqs;
		
	}

	public static Alignment importAlignment(String dataDir, String fileName) throws Exception{
		
		String filePath = dataDir+fileName;
		
		BufferedReader in = new BufferedReader(new FileReader(filePath));
		FastaImporter importer = new FastaImporter(in, Nucleotides.INSTANCE);
		Alignment alignment = importer.importAlignment();

		return alignment;
	}
	
	public static Alignment importShortReads(String dataDir, String fileName) throws Exception{
		
		String filePath = dataDir+fileName;
		
		BufferedReader in = new BufferedReader(new FileReader(filePath));
		ShortReadImporter importer = new ShortReadImporter(in, ShortReads.INSTANCE);
		Alignment alignment = importer.importAlignment();

		return alignment;
	}
	
		
	public static Tree importTree(String dataDir, String fileName) throws Exception{
		
		String filePath = dataDir+fileName;

		BufferedReader in = new BufferedReader(new FileReader(filePath));
		NewickImporter importer = new NewickImporter(in);
		Tree tree = importer.importTree(null);

		return tree;
	}
	

	
	
}
