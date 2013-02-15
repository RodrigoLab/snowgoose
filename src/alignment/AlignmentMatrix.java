package alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.CharSequenceUtils;
import org.apache.commons.lang3.CharSetUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.analysis.solvers.RegulaFalsiSolver;

import com.google.common.base.Strings;
import com.google.common.collect.ArrayTable;
import com.google.common.primitives.Chars;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;


/*
Represent the full alignment with a set of short reads


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
public class AlignmentMatrix {
	
	
	private static final char GAP = '-';
	private static Random rand = new Random();
	
	char[][] matrix;
//	ArrayList<char[]> matrix;
	int noHap;
	int length;
	AlignmentMapping aMap;
	
	SimpleAlignment alignment;
	private int[] oldhap;
	
//	public AlignmentMatrix(int noSeq, int length) {
//		this.noSeq = noSeq;
//		this.length = length;
////		matrix = new ArrayList<>();
////		matrix = new char[]
//		
//	}
	
	public AlignmentMatrix(AlignmentMapping aMap, int noHap) {
		this.aMap = aMap;
		this.length = aMap.getLength();
		this.noHap = noHap;
		matrix = new char[this.noHap][length];
		
		alignment = new SimpleAlignment();
		alignment.setDataType(Nucleotides.INSTANCE);
		init();
	}

	public AlignmentMatrix(AlignmentMapping aMap, Alignment trueAlignment) {
		this(aMap, trueAlignment.getSequenceCount());
		alignment = (SimpleAlignment) trueAlignment;
		for (int i = 0; i < noHap; i++) {
			Sequence s = trueAlignment.getSequence(i);
//			for (int j = 0; j < length; j++) {
				s.getChars(0, length, matrix[i], 0);
//			}
//			
		}
	
	}

	private void init(){
		for (int i = 0; i < noHap; i++) {
			randomSeq(i);
		}
	}
	public void randomSeq(int hapIndex) {
		randomSeq(hapIndex, 0, length);
		
	}
	
	
	public void randomSeq(int hapIndex, int start, int end){
		
		for (int p = start; p < end; p++) {
			int size = aMap.mapToSrp[p].size();
	
			if (size != 0){
				int srpIndex = aMap.mapToSrp[p].get(rand.nextInt(size));
//				matrix[hapIndex][p] = aMap.getShortRead(index).getFullSrpCharAt(p);
				swapBase(hapIndex, p, srpIndex);
			}
			else{
				matrix[hapIndex][p] = GAP;
			}
		}
	}

	public String getHaplotype(int i){
		
		String hap = String.valueOf(matrix[i]);
		return hap;
	}
	
	public void swapSrp(int hapIndex, int start, int end, int srpIndex){
		String srp = aMap.getFullSrp(srpIndex);
		for (int p = start; p < end; p++) {
			matrix[hapIndex][p] = srp.charAt(p);
		}
	}
	
	public void swapBase() {
		
		int hapIndex = rand.nextInt(noHap);
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

			int size = aMap.mapToSrp[pos].size();
			if (size != 0){
				int srpIndex = aMap.mapToSrp[pos].get(rand.nextInt(size));
	//			matrix[hapIndex][pos] = c;
				swapBase(hapIndex, pos, srpIndex);
			}
			else{
				matrix[hapIndex][pos] = GAP;
			}
		
		
	}
	private void swapBase(int hapIndex, int pos, int srpIndex){
		oldhap = new int[]{hapIndex, pos, matrix[hapIndex][pos]};
		char c = aMap.getShortReadCharAt(srpIndex, pos);
		
		matrix[hapIndex][pos] = c;
//		System.out.println(c +"\t"+ matrix[hapIndex][pos] +"\t"+ oldhap[2]);
		
	}

	//	public void addSeq(int index, char[] seqs) {
//		if(seqs.length != length){
//			System.err.println("Error, different length\t"+seqs.length +"\t"+ length);
//			System.exit(-1);
//		}
//		
//		for (int i = 0; i < length; i++) {
//			matrix[index][i] = seqs[i]; 
//		}
//		
//	}
	public String testGetSeq(){
		StringBuilder seq = new StringBuilder();
//		for (int i = 0; i < am.map.length; i++) {
		for (int i = 0; i <length; i++) {
			int size = aMap.mapToSrp[i].size();
//			System.out.print(i +"\t"+ am.mapToSrp[i]+"\t");
			if (size != 0){
				int index = aMap.mapToSrp[i].get(rand.nextInt(size));
				seq.append( aMap.getShortRead(index).getFullSrpCharAt(i) );
				
//				System.out.print(index+"\t"+ am.srp.get(index).charAt(i) +"\t");
//				for (int j = 0; j < am.map[i].size(); j++) {
//					System.out.print(am.map[i].get(j) +":"+am.srp.get(am.map[i].get(j)).charAt(i)+"\t" );
//				}
			}
			else{
				seq.append(GAP);
				
			}

		}
		System.out.println(seq.toString());
		return seq.toString();
	}
	
	public void testMultipleSeq(){
		
		for (int i = 0; i < matrix.length; i++) {
			
			randomSeq(i,0,length);
			System.out.println(String.valueOf(matrix[i]));
		}
		testSwapping();
	}
	
	public void testSwapping(){
		System.out.println("swapping");
		swapSrp(0,700,800,1);
		swapSrp(0,800,900,2);
		System.out.println(String.valueOf(matrix[0]));
	}
	
	public void toAlignment(){
		
		System.out.println("toAlignment");
		
		for (int i = 0; i < matrix.length; i++) {
			Sequence seq = new Sequence(String.valueOf(matrix[i]));
			
			alignment.addSequence(seq);	
//			System.out.println(seq.getDataType().getDescription());
			
		}
//		0
//		nucleotide
//		for (int i = 0; i < alignment.getSequenceCount(); i++) {
//			System.out
//					.println(alignment.getSequence(i).getSequenceString());
//		}
		int seqIndex = 1;
		System.out.println(alignment.getSequence(seqIndex).getSequenceString());
		

		Sequence s = alignment.getSequence(seqIndex);
		swapSrp(seqIndex ,800,1000,0);
		s.setSequenceString(String.valueOf(matrix[seqIndex]));
		//TODO setSequenceString is slowish
		System.out.println(alignment.getSequence(seqIndex).getSequenceString());
		
		
		Alignment aa = alignment;
		
	}

	public char[][] getCharMatrix() {
		return matrix;
	}

	public int getLength() {
		
		return length;
	}
	
	public int getNoHap(){
		return noHap;
	}

	public void reject() {

		matrix[oldhap[0]][oldhap[1]] = (char) oldhap[2];

	}

	
	
			
}