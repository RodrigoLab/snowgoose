package srp.spectrum;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import dr.evolution.alignment.Alignment;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;

public class SpectrumLogger extends MCLogger {

	

	private Alignment trueAlignment;
	


	public SpectrumLogger(SpectrumAlignmentModel alignment, Alignment trueAlignment, String fileName,
			int logEvery) throws IOException {
		this(alignment, fileName, logEvery);
		this.trueAlignment = trueAlignment;
		
	}
	

	protected SpectrumAlignmentModel spectrum;
	protected String[] taxaList;
	protected int count;


	public SpectrumLogger(SpectrumAlignmentModel alignment, String fileName, int logEvery) throws IOException {
		this(alignment, new TabDelimitedFormatter(new PrintWriter(new FileWriter(fileName))),logEvery);
//		TabDelimitedFormatter formatter = new TabDelimitedFormatter(new PrintWriter(new FileOutputStream(new File(fileName))));
//		TabDelimitedFormatter formatter = new TabDelimitedFormatter(new PrintWriter(new FileWriter(fileName)));
	}

	public SpectrumLogger(SpectrumAlignmentModel alignment, LogFormatter formatter, int logEvery) {
		super(formatter, logEvery, false);
		this.spectrum = alignment;
		count = alignment.getSpectrumCount();
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
				buffer.append(taxaList[i]).append("\n")
						.append(spectrum.getSpectrumString(i))
						.append("\n");
			}
			logLine(buffer.toString());
    	}

    }

    @Override
	public void stopLogging() {
    
        super.stopLogging();
    }
    
    //Contain true alignment? how is that useful? 
    //calculate the frequency delta?
	public void logTrue(long state) {

    	if (logEvery > 0 && (state % logEvery == 0)) {
//	    		(logEvery < 0 || ((state % logEvery) == 0));
    		
	        StringBuffer buffer = new StringBuffer("Alignment STATE_");
	        buffer.append(state).append("\n");
			for (int i = 0; i < count; i++) {
				buffer.append(taxaList[i])
//						.append(spectrum.getAlignedSequenceString(i))
						.append("\n");
			}
//			buffer.append(HaplotypeModelUtils.calculeteSPSArrayString(spectrum, trueAlignment));
			logLine(buffer.toString());
    	}

    
    }
}