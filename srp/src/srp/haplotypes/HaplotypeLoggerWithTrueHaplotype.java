package srp.haplotypes;

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
			for (int i = 0; i < count; i++) {
				buffer.append(taxaList[i])
						.append(alignment.getAlignedSequenceString(i))
						.append("\n");
			}
			buffer.append(HaplotypeModelUtils.calculeteSPSArrayString(alignment, trueAlignment));
			logLine(buffer.toString());
    	}

    
    }
}