package srp.evolution.haplotypes;

import java.io.IOException;

import dr.evolution.alignment.Alignment;

public class HaplotypeLoggerWithTrueHaplotype extends HaplotypeLogger{
	
	private Alignment trueAlignment;
	


	public HaplotypeLoggerWithTrueHaplotype(Alignment alignment, Alignment trueAlignment, String fileName,
			int logEvery) throws IOException {
		super(alignment, fileName, logEvery);
		this.trueAlignment = trueAlignment;
		
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
			logLine(buffer.toString());
    	}

    
    }
}