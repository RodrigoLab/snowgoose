package io;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.io.Importer;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.SequenceList;
import dr.evolution.sequence.Sequences;
import dr.evolution.util.Taxon;

import java.io.EOFException;
import java.io.IOException;
import java.io.Reader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.StringTokenizer;


public class ShortReadImporter  extends Importer {

public static final char FASTA_FIRST_CHAR = '>';
    

    

	
    /**
     * Constructor
     */
	public ShortReadImporter(Reader reader) {
		super(reader);
		// TODO Auto-generated constructor stub
	}

//    public FastaImporter(Reader reader, Writer commentWriter, DataType dataType) {
//        super(reader, commentWriter);
//        setCommentDelimiters('\0', '\0', '\0');
//
//        this.dataType = dataType;
//    }

    /**
     * importAlignment.
     */
    public Sequences importAlignment() throws IOException, ImportException
    {
        Sequences allReads = new Sequences();  

        try {
            // find fasta line start
            while (read() != FASTA_FIRST_CHAR) {
            }

            do {
                final String name = readLine().trim();
                StringBuffer seq = new StringBuffer();

                readSequence(seq, dataType, "" + FASTA_FIRST_CHAR, Integer.MAX_VALUE, "-", "?", "", "");

//                alignment.addSequence();
                allReads.addSequence(  new Sequence(new Taxon(name.toString()), seq.toString()) );

            } while (getLastDelimiter() == FASTA_FIRST_CHAR);
        } catch (EOFException e) {
            // catch end of file the ugly way.
        }

        return allReads;
    }

 
    private DataType dataType;
    private int maxNameLength = 10;

	
	
}
