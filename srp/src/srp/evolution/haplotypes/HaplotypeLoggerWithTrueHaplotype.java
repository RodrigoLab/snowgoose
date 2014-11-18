package srp.evolution.haplotypes;

import java.io.IOException;
import java.util.Iterator;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;

public class HaplotypeLoggerWithTrueHaplotype extends HaplotypeLogger{
	
	private Alignment trueAlignment;
	private SimpleAlignment doubleAlignment;


	public HaplotypeLoggerWithTrueHaplotype(Alignment alignment, Alignment trueAlignment, String fileName,
			int logEvery) throws IOException {
		super(alignment, fileName, logEvery);
		this.trueAlignment = trueAlignment;
//		this.doubleAlignment = (SimpleAlignment) trueAlignment;
		doubleAlignment = new SimpleAlignment();
		
		for (int i = 0; i < alignment.getSequenceCount(); i++) {
			doubleAlignment.addSequence(alignment.getSequence(i));
		}
		for (int i = 0; i < trueAlignment.getSequenceCount(); i++) {
			doubleAlignment.addSequence(trueAlignment.getSequence(i));
		}
//		System.out.println(doubleAlignment.getSequenceCount());
		
	}
	@Override
	public void log(long state) {
		
    	if (logEvery > 0 && (state % logEvery == 0)) {
//	    		(logEvery < 0 || ((state % logEvery) == 0));
    		
	        StringBuffer buffer = new StringBuffer("Alignment STATE_");
	        buffer.append(state).append("\n");
	        if(state == 0){
	        	buffer.append(SPSDist.calculeteSPSArrayString(trueAlignment, trueAlignment));
				
			}
			for (int i = 0; i < count; i++) {
				buffer.append(taxaList[i])
						.append(alignment.getAlignedSequenceString(i))
						.append("\n");
			}
			buffer.append(SPSDist.calculeteSPSArrayString(alignment, trueAlignment));
			buffer.append("\n");
			buffer.append(SPSDist.calculeteSPSArrayString(alignment, alignment));
			buffer.append("[DoubleMatrix]\n");
			buffer.append(SPSDist.calculeteSPSArrayStringFormat(doubleAlignment, doubleAlignment));
			logLine(buffer.toString());
    	}

    
    }
}