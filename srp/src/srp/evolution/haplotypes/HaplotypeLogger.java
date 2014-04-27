package srp.evolution.haplotypes;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import dr.evolution.alignment.Alignment;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;

public class HaplotypeLogger extends MCLogger {
//
//	public HaplotypeLogger(LogFormatter formatter, int logEvery,
//			boolean performanceReport, int performanceReportDelay) {
//		super(formatter, logEvery, performanceReport, performanceReportDelay);
//	}
	
	protected Alignment alignment;
	protected String[] taxaList;
	protected int count;


	public HaplotypeLogger(Alignment alignment, String fileName, int logEvery) throws IOException {
		this(alignment, new TabDelimitedFormatter(new PrintWriter(new FileWriter(fileName))),logEvery);
//		TabDelimitedFormatter formatter = new TabDelimitedFormatter(new PrintWriter(new FileOutputStream(new File(fileName))));
//		TabDelimitedFormatter formatter = new TabDelimitedFormatter(new PrintWriter(new FileWriter(fileName)));
	}

	public HaplotypeLogger(Alignment alignment, LogFormatter formatter, int logEvery) {
		super(formatter, logEvery, false);
		this.alignment = alignment;
		count = alignment.getSequenceCount();
		taxaList = new String[count];
		for (int i = 0; i < taxaList.length; i++) {
			taxaList[i] = ">"+alignment.getTaxonId(i)+"\t";
		}
	}


    @Override
	public void startLogging() {
    	logLine("START!!");
    }

    @Override
	public void log(long state) {

    	if (logEvery > 0 && (state % logEvery == 0)) {
//    		(logEvery < 0 || ((state % logEvery) == 0));
    		
	        StringBuffer buffer = new StringBuffer("Alignment STATE_");
	        buffer.append(state).append("\n");
			for (int i = 0; i < count; i++) {
				buffer.append(taxaList[i])
						.append(alignment.getAlignedSequenceString(i))
						.append("\n");
			}
			logLine(buffer.toString());
    	}

    }

    @Override
	public void stopLogging() {
    
        super.stopLogging();
    }

}