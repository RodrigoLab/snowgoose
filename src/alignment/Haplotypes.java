package alignment;

import java.util.Arrays;
import java.util.Random;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import dr.util.NumberFormatter;


/*
Represent the full alignment with a set of short reads

http://www.ncbi.nlm.nih.gov/pmc/articles/PMC148477/pdf/272682.pdf
Alignment scores
 - sum-of-pairs score, might need mapping
 - column score, 1 or 0 column match

// resusable buffer for case 2

    
    char[] cbuff = new char[1024 * 1024];

    // a reusable field we will use for reflection in case 4
    Field field = String.class.getDeclaredField("value");
    field.setAccessible(true);


        // CASE 2 - Copy chars to a buffer
        time = System.currentTimeMillis();

        // check chars of string 10M times using cbuff[i] method
        for (int n = 0; n < 1_000_000; n++) {
            int count = data.length();
            data.getChars(0, count, cbuff, 0);
            for (int i = 0; i < count; i++) {
                if (cbuff[i] <= ' ') {
                    throw new IllegalDataException("Found whitespace");
                }
            }
        }



        // CASE 4 - Use reflection to access String char[]
        time = System.currentTimeMillis();

        // check chars of string 10M times using reflection method
        for (int n = 0; n < 1_000_000; n++) {
            final char[] chars = (char[]) field.get(data);
            final int len = chars.length;
            for (int i = 0; i < len; i++) {
                if (chars[i] <= ' ') {
                    throw new Exception("Found whitespace");
                }
            }
        }


*/
//public class Haplotypes implements Alignment{
public class Haplotypes {
	
	public static final int[] NULL_SWAPINFO = new int[4];
	public static final char GAP = '-';
	
	private static Random rand = new Random();
	
	char[][] matrix;
//	ArrayList<char[]> matrix;
	int haplotypesCount;
	int haplotypesLength;
	AlignmentMapping aMap;
	
	SimpleAlignment alignment;
	Taxon[] taxons;
	private int[] swapInfo = new int[NULL_SWAPINFO.length]; //TODO: Expand to a class later, use this to speed up likelihood calculation
	
//	public AlignmentMatrix(int noSeq, int length) {
//		this.noSeq = noSeq;
//		this.length = length;
////		matrix = new ArrayList<>();
////		matrix = new char[]
//		
//	}
	
	public Haplotypes(AlignmentMapping aMap, int noHap) {
		this.aMap = aMap;
		this.haplotypesLength = aMap.getLength();
		this.haplotypesCount = noHap;
		matrix = new char[this.haplotypesCount][haplotypesLength];
		
		alignment = new SimpleAlignment();
		alignment.setDataType(Nucleotides.INSTANCE);
		
		taxons = new Taxon[haplotypesCount];
		for (int i = 0; i < taxons.length; i++) {
			taxons[i] = new Taxon("hap"+(i+1)); //TODO get taxons name from tree, or at lesat sync the name
		}
		initSeqs();
		swapInfo = new int[4];
		
	}

	public Haplotypes(AlignmentMapping aMap, Alignment trueAlignment) {
		this(aMap, trueAlignment.getSequenceCount());
		alignment = (SimpleAlignment) trueAlignment;
		alignmentToMatrix();
	
	}
//	public void setAlignment(Alignment trueAlignment){
//		
//		alignment = (SimpleAlignment) trueAlignment;
//		
//	}
	private void alignmentToMatrix(){
		for (int i = 0; i < haplotypesCount; i++) {
			Sequence tempSeq = alignment.getSequence(i);
				tempSeq.getChars(0, haplotypesLength, matrix[i], 0);
		}
		
	}
	private void matrixToAlignment(){

//		System.out.println("toAlignment");
//		Sequence tempSeq = new Sequence();
		for (int i = 0; i < haplotypesCount; i++) {
			
//			alignment.addSequence(tempSeq);	
			Sequence tempSeq = alignment.getSequence(i);
			tempSeq.setSequenceString( getHaplotype(i) );
		}
//		0
//		nucleotide
//		for (int i = 0; i < alignment.getSequenceCount(); i++) {
//			System.out
//					.println(alignment.getSequence(i).getSequenceString());
////		}
//		int seqIndex = 1;
//		System.out.println(alignment.getSequence(seqIndex).getSequenceString());
//		
//
//		Sequence s = alignment.getSequence(seqIndex);
//		swapSrp(seqIndex ,800,1000,0);
//		s.setSequenceString(String.valueOf(matrix[seqIndex]));
//		//TODO setSequenceString is slowish
//		System.out.println(alignment.getSequence(seqIndex).getSequenceString());
//		
		
//		Alignment aa = alignment;
		
	}
	private void initSeqs(){
		for (int i = 0; i < haplotypesCount; i++) {
			randomSeq(i);
			Sequence s = new Sequence(taxons[i], getHaplotype(i));
			alignment.addSequence(s);
		}
	}
	
	public void randomSeq(int hapIndex) {
		randomSeq(hapIndex, 0, haplotypesLength);
		
	}
	
	
	public void randomSeq(int hapIndex, int start, int end){
		
		for (int p = start; p < end; p++) {
			swapBase(hapIndex, p);
//			int size = aMap.mapToSrp[p].size();
//	
//			if (size != 0){
//				int srpIndex = aMap.mapToSrp[p].get(rand.nextInt(size));
////				matrix[hapIndex][p] = aMap.getShortRead(index).getFullSrpCharAt(p);
//				swapBase(hapIndex, p, srpIndex);
//			}
//			else{
//				matrix[hapIndex][p] = GAP;
//			}
		}
	}

	public String getHaplotype(int i){

		return String.valueOf(matrix[i]);
	}
	
	public void swapSrp(int hapIndex, int start, int end, int srpIndex){
		String srp = aMap.getFullSrp(srpIndex);
		for (int p = start; p < end; p++) {
			matrix[hapIndex][p] = srp.charAt(p);
		}
	}
	
	public void swapBase() {
		
		int hapIndex = rand.nextInt(haplotypesCount);
		swapBase(hapIndex);
	}

	public void swapBase(int hapIndex){
		int pos = rand.nextInt(aMap.getLength());
		swapBase(hapIndex, pos);
//		int size = aMap.mapToSrp[pos].size();
//		int srpIndex = aMap.mapToSrp[pos].get(rand.nextInt(size));
////		char c = aMap.getShortReadCharAt(srpIndex, pos);
//		
//		swapBase(hapIndex, pos, srpIndex);
		
	}
	
	public void swapBase(int hapIndex, int pos){

		swapInfo[0] = hapIndex;
		swapInfo[1] = pos;
		swapInfo[2] = matrix[hapIndex][pos];
		
		int size = aMap.mapToSrp[pos].size();
		if (size != 0){
			int srpIndex = aMap.mapToSrp[pos].get(rand.nextInt(size));
			swapBase(hapIndex, pos, srpIndex);
		}
		else{
			matrix[hapIndex][pos] = GAP;
			swapInfo[3] = GAP;
		}
	
	
	}

	private void swapBase(int hapIndex, int pos, int srpIndex){

		char c = aMap.getShortReadCharAt(srpIndex, pos);
		matrix[hapIndex][pos] = c;
		swapInfo[3] = c;

		// System.out.println(c +"\t"+ matrix[hapIndex][pos] +"\t"+ oldhap[2]);
		
	}
	private void swapBase(int hapIndex, int pos, int srpIndex, boolean t){
//		swapInfo = new int[]{hapIndex, pos, matrix[hapIndex][pos]};
		char c = 0;//matrix[hapIndex][pos];
		do{
			c = aMap.getShortReadCharAt(srpIndex, pos);
//			System.out.println("A" +"\t"+ matrix[hapIndex][pos] +"\t"+  c);
		} while (matrix[hapIndex][pos] == c );
		matrix[hapIndex][pos] = c;
		swapInfo[3] = c;
//		System.out.println(c +"\t"+ matrix[hapIndex][pos] +"\t"+ oldhap[2]);
		
	}

	public double calculateSPS(){
		double sps = 0;
		for (int i = 0; i < haplotypesCount; i++) {
//			System.out.println(String.valueOf(matrix[i]));
			for (int j = 0; j < i; j++) {
//				System.out.println(String.valueOf(matrix[j]));
				for (int b = 0; b < haplotypesLength; b++) {
					int c = matrix[i][b] - matrix[j][b];
					sps += (c==0)? 0: 1;
//					if(c!=0){
//						sps+=1;
//						System.out.println(i +"\t"+ j +"\t"+ b +"\t"+ matrix[i][b] +"\t"+  matrix[j][b] +"\t"+ sps);
//					}
				}
				
			}
		}
		return sps;
	}


	public char[][] getCharMatrix() {
		return matrix;
	}

	public int getLength() {
		
		return haplotypesLength;
	}
	
	public int getHaplotypesCount(){
		return haplotypesCount;
	}

	public void reject() {

		matrix[swapInfo[0]][swapInfo[1]] = (char) swapInfo[2];

	}

	
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
        for (int i = 0; i < getSequenceCount(); i++) {
            String name = formatter.formatToFieldWidth(getTaxonId(i), 10);
            buffer.append(">" + name + "\n");
            buffer.append(getAlignedSequenceString(i) + "\n");
        }

        return buffer.toString();
    
		
	}

//	@Override
	public int getSequenceCount() {
		return getHaplotypesCount();
	}

//	@Override
//	public Sequence getSequence(int i) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public void setSequenceAttribute(int index, String name, Object value) {
//		// TODO Auto-generated method stub
//		
//	}
//
//	@Override
//	public Object getSequenceAttribute(int index, String name) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public int getTaxonCount() {
//		// TODO Auto-generated method stub
//		return 0;
//	}
//
//	@Override
//	public Taxon getTaxon(int taxonIndex) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
	public String getTaxonId(int taxonIndex) {

		return taxons[taxonIndex].getId();
	}
//
//	@Override
//	public int getTaxonIndex(String id) {
//		// TODO Auto-generated method stub
//		return 0;
//	}
//
//	@Override
//	public int getTaxonIndex(Taxon taxon) {
//		// TODO Auto-generated method stub
//		return 0;
//	}
//
//	@Override
//	public List<Taxon> asList() {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public Object getTaxonAttribute(int taxonIndex, String name) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public String getId() {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public void setId(String id) {
//		// TODO Auto-generated method stub
//		
//	}
//
//	@Override
//	public Iterator<Taxon> iterator() {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public int getSiteCount() {
//		// TODO Auto-generated method stub
//		return 0;
//	}
//
//	@Override
//	public int[] getSitePattern(int siteIndex) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public int getPatternIndex(int siteIndex) {
//		// TODO Auto-generated method stub
//		return 0;
//	}
//
//	@Override
//	public int getState(int taxonIndex, int siteIndex) {
//		// TODO Auto-generated method stub
//		return 0;
//	}
//
//	@Override
//	public int getPatternCount() {
//		// TODO Auto-generated method stub
//		return 0;
//	}
//
//	@Override
//	public int getStateCount() {
//		// TODO Auto-generated method stub
//		return 0;
//	}
//
//	@Override
//	public int getPatternLength() {
//		// TODO Auto-generated method stub
//		return 0;
//	}
//
//	@Override
//	public int[] getPattern(int patternIndex) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public int getPatternState(int taxonIndex, int patternIndex) {
//		// TODO Auto-generated method stub
//		return 0;
//	}
//
//	@Override
//	public double getPatternWeight(int patternIndex) {
//		// TODO Auto-generated method stub
//		return 0;
//	}
//
//	@Override
//	public double[] getPatternWeights() {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public DataType getDataType() {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public double[] getStateFrequencies() {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public void setDataType(DataType dataType) {
//		// TODO Auto-generated method stub
//		
//	}
//
//	@Override
	public String getAlignedSequenceString(int sequenceIndex) {
		
		return getHaplotype(sequenceIndex);
	}
//
//	@Override
//	public String getUnalignedSequenceString(int sequenceIndex) {
//		// TODO Auto-generated method stub
//		return null;
//	}

	

	public static int[][] calculeteSPSArray(Haplotypes h1, Haplotypes h2){
		
		
		char[][] m1 = h1.getCharMatrix();
		char[][] m2 = h2.getCharMatrix();
		int hapLength = m1.length;
		int seqLength = m1[0].length;
		int sps[][] = new int[hapLength][hapLength];
		if (seqLength != m2[0].length){
			System.err.println("Incompariable alignments lenght: "+m1[0].length +" and "+  m2[0].length);
		}
		for (int i = 0; i < m1.length; i++) {
			for (int j = 0; j < m2.length; j++) {
				sps[i][j] = 0;
				for (int l = 0; l < seqLength; l++) {
//					System.out.println(((m1[i][l] - m2[j][l]) == 0) +"\t"+ m1[i][l] +"\t"+  m2[j][l]);
//					sps +=  ((m1[i][l] - m2[j][l]) == 0)? 0:1;
					sps[i][j] +=  ((m1[i][l] - m2[j][l]) == 0)? 0:1;
					
				}
				
			}
		}
		return sps;
		
	}
	public static int calculeteSPS(Haplotypes h1, Haplotypes h2){
		
		int sps = 0;
		char[][] m1 = h1.getCharMatrix();
		char[][] m2 = h2.getCharMatrix();
		int seqLength = m1[0].length;
		if (seqLength != m2[0].length){
			System.err.println("Incompariable alignments lenght: "+m1[0].length +" and "+  m2[0].length);
		}
		for (int i = 0; i < m1.length; i++) {
			for (int j = 0; j < m2.length; j++) {

				for (int l = 0; l < seqLength; l++) {
//					System.out.println(((m1[i][l] - m2[j][l]) == 0) +"\t"+ m1[i][l] +"\t"+  m2[j][l]);
//					sps +=  ((m1[i][l] - m2[j][l]) == 0)? 0:1;
					sps +=  ((m1[i][l] - m2[j][l]) == 0)? 0:1;
					
				}
				
			}
		}
		return sps;
		
	}

	public int[] getSwapInfo() {
		return swapInfo;
	}

	public Alignment getAlignment() {
		matrixToAlignment();
		return alignment;
	}
	
}