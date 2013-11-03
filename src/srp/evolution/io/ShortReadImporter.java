package srp.evolution.io;

import java.io.EOFException;
import java.io.IOException;
import java.io.Reader;
import java.io.Writer;

import srp.evolution.datatype.ShortReads;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.io.FastaImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

public class ShortReadImporter extends FastaImporter {

	private final static DataType DATA_TYPE= ShortReads.INSTANCE;
	
	public ShortReadImporter(Reader reader) {
		super(reader, DATA_TYPE);
	}
	public ShortReadImporter(Reader reader, DataType dataType) {
		super(reader, DATA_TYPE);
		// TODO Auto-generated constructor stub
	}

	public ShortReadImporter(Reader reader, Writer commentWriter,
			DataType dataType) {
		super(reader, commentWriter, dataType);
		// TODO Auto-generated constructor stub
	}

    public Alignment importAlignment() throws IOException, ImportException
    {
        SimpleAlignment alignment = null;

        try {
            // find fasta line start
            while (read() != FASTA_FIRST_CHAR) {
            }

            do {
                final String name = readLine().trim();
                StringBuffer seq = new StringBuffer();

                readSequence(seq, null, "" + FASTA_FIRST_CHAR, Integer.MAX_VALUE, "-", "?", "", "");

                if (alignment == null) {
                    alignment = new SimpleAlignment();
                    alignment.setDataType(DATA_TYPE);
                }
                Sequence sequence = new Sequence(new Taxon(name.toString()), seq.toString());
                sequence.setDataType(ShortReads.INSTANCE);
                alignment.addSequence(sequence);

            } while (getLastDelimiter() == FASTA_FIRST_CHAR);
        } catch (EOFException e) {
            // catch end of file the ugly way.
        }

        return alignment;
    }

}
