package srp.evolution.haplotypes;

import java.util.Arrays;

import dr.evolution.alignment.Alignment;
import dr.util.NumberFormatter;

public class SPSDist {

	public static int calculeteSPS(Alignment h1, Alignment h2){
	//		
			int sps = 0;
			int[][] spsArray = calculateSPSCore(h1, h2);
			
			for (int i = 0; i < spsArray.length; i++) {
				for (int j = 0; j < spsArray[i].length; j++) {
					sps += spsArray[i][j];
				}
			}
			return sps;
			
		}

	public static int[][] calculeteSPSArray(Alignment h1, Alignment h2){
		
		int[][] sps = calculateSPSCore(h1,h2);
		return sps;
		
	}
//	return formatter.formatToFieldWidth(formatter.formatDecimal(acceptanceProb, 6), 11) + " ";//TODO SW Change
	protected static NumberFormatter formatter = new NumberFormatter(8);
	
	public static String calculeteSPSArrayStringFormat(Alignment h1, Alignment h2){
		int[][] sps = calculeteSPSArray(h1,h2);
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < sps.length; i++) {
			for (int j = 0; j < sps[i].length; j++) {
				
			sb.append( formatter.formatToFieldWidth(""+sps[i][j], 4) )
				.append(" ");
			
//			sb.append(Arrays.toString(sps[i])).append("\n");
			}
			sb.append("\n");
		}
		return sb.toString();
	}

	public static String calculeteSPSArrayString(Alignment h1, Alignment h2){
		int[][] sps = calculeteSPSArray(h1,h2);
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < sps.length; i++) {
			sb.append(Arrays.toString(sps[i])).append("\n");
		}
		return sb.toString();
	}

	private static int caluclateSPSSingle(String s1, String s2){
		int sps = 0;
		if (s1.length() != s2.length()){
			System.err.println("Incompariable sequence lenght: "+ s1.length() +" and "+  s2.length());
		}
		for (int i = 0; i < s1.length(); i++) {
			int c = s1.charAt(i) - s2.charAt(i);
			sps += (c == 0)? 0:1;
		}
	
		return sps;
	}

	private static int[][] calculateSPSCore(Alignment h1, Alignment h2) {
		int sps[][];
		int hapCount = h1.getSequenceCount();
		int hapCount2 = h2.getSequenceCount();
		int seqLength = h1.getSiteCount();

		if (seqLength != h2.getSiteCount()) {
			System.err.println("Incompariable alignments lenght: " + seqLength
					+ " and " + h2.getSiteCount());
			System.exit(-1);
		}
		
		if (hapCount != hapCount2) {
			
			String[] s1 = new String[hapCount];
			String[] s2 = new String[hapCount2];
			
			for (int i = 0; i < s1.length; i++) {
				s1[i] = h1.getAlignedSequenceString(i);
			}
			for (int i = 0; i < s2.length; i++) {
				s2[i] = h2.getAlignedSequenceString(i);
			}
			
			sps = new int[hapCount][hapCount2];
			for (int i = 0; i < s1.length; i++) {
				for (int j = 0; j < s2.length; j++) {
					sps[i][j] = caluclateSPSSingle(s1[i], s2[j]);
				}
			}

		} else {
			sps = new int[hapCount][hapCount];
			String[] s1 = new String[hapCount];
			String[] s2 = new String[hapCount];
			for (int i = 0; i < hapCount; i++) {
				s1[i] = h1.getAlignedSequenceString(i);
				s2[i] = h2.getAlignedSequenceString(i);
			}
			for (int i = 0; i < hapCount; i++) {
				for (int j = 0; j < hapCount; j++) {
					sps[i][j] = caluclateSPSSingle(s1[i], s2[j]);
				}
			}

		}
		return sps;
	}


}
