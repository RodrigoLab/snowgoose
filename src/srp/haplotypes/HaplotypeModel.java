package srp.haplotypes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import dr.util.Attributable;
import dr.util.NumberFormatter;


/*


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
public class HaplotypeModel implements Alignment{
//public class HaplotypeModel extends Abstract{
	
	public static final int[] NULL_SWAPINFO = new int[4];
	public static final char GAP = '-';
	public static final String TAXON_PREFIX = "Hap_";
	
	
	private static Random rand = new Random();


    private Attributable.AttributeHelper attributes = null;

	private DataType dataType = null;
	
	char[][] matrix;

//	int haplotypesCount;
	int haplotypesLength;
	AlignmentMapping aMap;
	
	ArrayList<Haplotype> haplotypes;
//	SimpleAlignment alignment;
//	Taxa taxa;	//FIXME might be better/easier to like with trees, or might not need at all
//	Taxon[] taxons;
	private int[] swapInfo = new int[NULL_SWAPINFO.length]; 
	//TODO: Expand to a class later, use this to speed up likelihood calculation
	
	
	
	
	
	
	
	
	
	
	private void setupData(AlignmentMapping aMap) {
		this.aMap = aMap;
		this.haplotypesLength = this.aMap.getLength();
		
//		this.haplotypesCount = hapCount;

//		taxa = new Taxa();
		
//		taxons = new Taxon[haplotypesCount];
//		alignment = new SimpleAlignment();
		
		haplotypes = new ArrayList<Haplotype>();

		setDataType(Nucleotides.INSTANCE);
		
	}
	
	private void setupAlignment(Alignment trueAlignment) {

		for (int i = 0; i < trueAlignment.getSequenceCount(); i++) {
			Haplotype haplotype = new Haplotype(trueAlignment.getSequence(i));
			addHaplotype(haplotype);
		}
		matrix = new char[getHaplotypeCount()][getHaplotypeLength()];
		
    }
	
	private void initSeqs(int hapCount){
		matrix = new char[hapCount][getHaplotypeLength()];
		char[] temp = new char[getHaplotypeLength()];
		Arrays.fill(temp, DataType.GAP_CHARACTER);
		String tempSeq = String.valueOf(temp);
		
		for (int i = 0; i < hapCount; i++) {
			Taxon t = new Taxon(TAXON_PREFIX+i); 
			Haplotype haplotype = new Haplotype(t, tempSeq);
			addHaplotype(haplotype);

//			String tempSeq = 
					randomSeq(i);

		}
	}
	private void addHaplotype(Haplotype haplotype) {
		haplotype.setDataType(dataType);
	    haplotypes.add(haplotype);
	
	}

	public HaplotypeModel(AlignmentMapping aMap, int hapCount) {
		setupData(aMap);
				
		initSeqs(hapCount);
		swapInfo = new int[4];//TODO remove later
		
	}

	public HaplotypeModel(AlignmentMapping aMap, Alignment trueAlignment) {
		setupData(aMap);
		setupAlignment(trueAlignment);
		alignmentToMatrix();
	
	}
	@Deprecated
	private void alignmentToMatrix(){

		for (int i = 0; i < getHaplotypeCount(); i++) {
			Haplotype tempSeq = haplotypes.get(i);
			tempSeq.getChars(0, getHaplotypeLength(), matrix[i], 0);
		}
		
	}
	
	public String randomSeq(int hapIndex) {
		String tempSeq = randomSeq(hapIndex, 0, getHaplotypeLength());
		return tempSeq;
	}
	
	
	public String randomSeq(int hapIndex, int start, int end){
		
		for (int p = start; p < end; p++) {
			swapBase(hapIndex, p);

		}
		String tempSeq = String.valueOf(matrix[hapIndex]);
		return tempSeq;
		
	}

	public void swapBase() {
		
		int hapIndex = rand.nextInt(getHaplotypeCount());
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
		
		char c = GAP;
		int size = aMap.mapToSrp[pos].size();
		if (size != 0) {
			int srpIndex = aMap.mapToSrp[pos].get(rand.nextInt(size));
			c = aMap.getShortReadCharAt(srpIndex, pos);
		}
		swapBase(hapIndex, pos, c);
	
	
	}

	private void swapBase(int hapIndex, int pos, char c){
		swapInfo[0] = hapIndex;
		swapInfo[1] = pos;
		swapInfo[2] = matrix[hapIndex][pos];
		swapInfo[3] = c;
		
		matrix[hapIndex][pos] = c;

		getHaplotype(hapIndex).setCharAt(pos, c);

		
	}



	public void swapSrp(int hapIndex, int start, int end, int srpIndex){
		String srp = aMap.getSrpFull(srpIndex);
		for (int p = start; p < end; p++) {
			matrix[hapIndex][p] = srp.charAt(p);
		}
	}

	public char[][] getCharMatrix() {
		return matrix;
	}

	public void reject() {

		matrix[swapInfo[0]][swapInfo[1]] = (char) swapInfo[2];
		getHaplotype(swapInfo[0]).setCharAt(swapInfo[1], (char)swapInfo[2]);
	}

	@Deprecated
	public Alignment getAlignment() {
		SimpleAlignment alignment = new SimpleAlignment();
		for (int i = 0; i < getHaplotypeCount(); i++) {
			Haplotype tempSeq = haplotypes.get(i);
			tempSeq.setHaplotypeString( getHaplotypeString(i) );
			alignment.addSequence(tempSeq);
		}
		
		return alignment;
	}

	
	public int calculateSPS(){
		int sps = 0;
		for (int i = 0; i < getHaplotypeCount(); i++) {
			for (int j = 0; j < i; j++) {
				for (int b = 0; b < getHaplotypeLength(); b++) {
					int c = matrix[i][b] - matrix[j][b];
					sps += (c==0)? 0: 1;
				}
			}
		}
		return sps;
	}

	public int[] getSwapInfo() {
		return swapInfo;
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
        for (int i = 0; i < getHaplotypeCount(); i++) {
            String name = formatter.formatToFieldWidth(getTaxon(i).getId(), 10);
            buffer.append(">" + name + "\n");
            buffer.append(getAlignedSequenceString(i) + "\n");
        }

        return buffer.toString();
    
		
	}


	public String getHaplotypeString(int i) {
	
//		return String.valueOf(matrix[i]);
		return getHaplotype(i).getSequenceString();
		
	}

	public int getHaplotypeLength() {
	
		return haplotypesLength;
	}

	public int getHaplotypeCount() {
		// TODO OR sequence.length()? change number of haplotypes
		return haplotypes.size();
	}

	public Haplotype getHaplotype(int i){
		return haplotypes.get(i);
	}

    // **************************************************************
    // SequenceList IMPLEMENTATION
    // **************************************************************
	
	/*
	 * Call getHaplotypesCount(), return haplotypesCount
	 */
	@Override
	public int getSequenceCount() {
		return getHaplotypeCount();
	}

	@Override
	public Sequence getSequence(int i) {
		// TODO think again
		return getHaplotype(i);
	}

	
	public void setSequenceAttribute(int index, String name, Object value) {
		// XXX double check
		Sequence sequence = getSequence(index);
        sequence.setAttribute(name, value);
	}

	
	public Object getSequenceAttribute(int index, String name) {
		//XXX double check
		Sequence sequence = getSequence(index);
        return sequence.getAttribute(name);
	}

    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

	@Override
	public int getTaxonCount() {
		return getHaplotypeCount();
	}
	
	@Override
	public Taxon getTaxon(int taxonIndex) {
		return getHaplotype(taxonIndex).getTaxon();
	}

	/**
     * @return the ID of the taxon of the ith sequence. If it doesn't have
     *         a taxon, returns the ID of the sequence itself.
     */
    public String getTaxonId(int taxonIndex) {
		//XXX double check
        Taxon taxon = getTaxon(taxonIndex);
        if (taxon != null)
            return taxon.getId();
        else
            throw new IllegalArgumentException("Illegal taxon index:" + taxonIndex);
    }

    /**
     * returns the index of the taxon with the given id.
     */
    @Override
	public int getTaxonIndex(String id) {
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            if (getTaxonId(i).equals(id)) return i;
        }
        return -1;
    }

    /**
     * returns the index of the given taxon.
     * must be the same object
     */
    @Override
	public int getTaxonIndex(Taxon taxon) {
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            if (getTaxon(i) == taxon) return i;
        }
        return -1;
    }

    @Override
	public List<Taxon> asList() {
        List<Taxon> taxa = new ArrayList<Taxon>();
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            taxa.add(getTaxon(i));
        }
        return taxa;
    }

	
//    /**
//     * Sets an named attribute for the taxon of a given sequence. If the sequence
//     * doesn't have a taxon then the attribute is added to the sequence itself.
//     *
//     * @param taxonIndex the index of the taxon whose attribute is being set.
//     * @param name       the name of the attribute.
//     * @param value      the new value of the attribute.
//     */
//    public void setTaxonAttribute(int taxonIndex, String name, Object value) {
//        Taxon taxon = getTaxon(taxonIndex);
//        if (taxon != null)
//            taxon.setAttribute(name, value);
//        else
//            setSequenceAttribute(taxonIndex, name, value);
//    }

    /**
     * @param taxonIndex the index of the taxon whose attribute is being fetched.
     * @param name       the name of the attribute of interest.
     * @return an object representing the named attributed for the given taxon.
     */
    @Override
	public Object getTaxonAttribute(int taxonIndex, String name) {
    	Taxon taxon = getTaxon(taxonIndex);
        if (taxon != null)
            return taxon.getAttribute(name);
        else
        	//TODO test this, might not needed
            return getSequenceAttribute(taxonIndex, name);
    }

    @Override
	public Iterator<Taxon> iterator() {
        return new Iterator<Taxon>() {
            private int index = -1;

            @Override
			public boolean hasNext() {
                return index < getTaxonCount() - 1;
            }

            @Override
			public Taxon next() {
                index++;
                return getTaxon(index);
            }

            @Override
			public void remove() { /* do nothing */ }
        };
    }


    // **************************************************************
    // SiteList IMPLEMENTATION
    // **************************************************************

    /**
     * Gets the pattern index at a particular site
     * @param siteIndex
     * @return siteIndex, identical to @param
     */
    @Override
	public int getPatternIndex(int siteIndex) {
        return siteIndex;
    }

    /**
     * @return number of sites - getHaplotypeLength()
     */
    @Override
	public int getSiteCount() {
        return getHaplotypeLength();
    }

    /**
     * Gets the pattern of site as an array of state numbers (one per sequence)
     *
     * @return the site pattern at siteIndex
     */
    @Override
	public int[] getSitePattern(int siteIndex) {

    	int n = getHaplotypeCount();

    	int[] pattern = new int[n];
        for (int i = 0; i < n; i++) {
            Haplotype hap = getHaplotype(i);

            if (siteIndex >= hap.getLength())
                pattern[i] = dataType.getGapState();
            else
                pattern[i] = hap.getState(siteIndex);
        }

        return pattern;
    }


    /**
     * @return the sequence state at (taxon, site)
     */
    @Override
	public int getState(int taxonIndex, int siteIndex) {
        Haplotype hap = getHaplotype(taxonIndex);

        if (siteIndex >= hap.getLength()) {
            return dataType.getGapState();
        }

        return hap.getState(siteIndex);
    }

//    /**
//     */
//    public void setState(int taxonIndex, int siteIndex, int state) {
//
//        Sequence seq = getSequence(taxonIndex);
//
//        if (siteIndex >= seq.getLength()) {
//            throw new IllegalArgumentException();
//        }
//
//        seq.setState(siteIndex, state);
//    }

    // **************************************************************
    // PatternList IMPLEMENTATION
    // **************************************************************

    /**
     * @return number of patterns
     */
    @Override
	public int getPatternCount() {
        return getSiteCount();
    }
 

    /**
     * @return number of states for this siteList
     */
    @Override
	public int getStateCount() {
        return getDataType().getStateCount();
    }

    /**
     * Gets the length of the pattern strings which will usually be the
     * same as the number of taxa
     *
     * @return the length of patterns
     */
    @Override
	public int getPatternLength() {
        return getSequenceCount();
    }

    /**
     * Gets the pattern as an array of state numbers (one per sequence)
     *
     * @return the pattern at patternIndex
     */
    public int[] getPattern(int patternIndex) {
    	//XXX double check
        return getSitePattern(patternIndex);
    }

    /**
     * @return state at (taxonIndex, patternIndex)
     */
    public int getPatternState(int taxonIndex, int patternIndex) {
    	//XXX double check
    	return getState(taxonIndex, patternIndex);
    }

    /**
     * Gets the weight of a site pattern (always 1.0)
     */
    @Override
	public double getPatternWeight(int patternIndex) {
        return 1.0;
    }

    /**
     * @return the array of pattern weights
     */
    @Override
	public double[] getPatternWeights() {
        double[] weights = new double[getSiteCount()];
        for (int i = 0; i < weights.length; i++)
            weights[i] = 1.0;
        return weights;
    }

	
    /**
     * @return the frequency of each state
     */
    public double[] getStateFrequencies() {
    	// XXX double check
        return PatternList.Utils.empiricalStateFrequencies(this);
    }

	@Override
	public DataType getDataType() {
		// XXX double check
		return dataType;
	}
	
    // **************************************************************
    // Alignment IMPLEMENTATION
    // **************************************************************

	@Override
	public void setDataType(DataType dataType) {
		// XXX double check
		this.dataType = dataType;
	}



	
	/**
	 * call getHaplotypeString(sequenceIndex);
	 */
	@Override
	public String getAlignedSequenceString(int sequenceIndex) {
		
		return getHaplotypeString(sequenceIndex);
	}

	/**
	 * call getHaplotypeString(sequenceIndex);
	 */
	@Override
	public String getUnalignedSequenceString(int sequenceIndex) {
		
		return getHaplotypeString(sequenceIndex);
	}
	
    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

	protected String id = null;

	/**
	 * @return the id.
	 */
	@Override
	public String getId() {
		return id;
	}

	/**
	 * Sets the id.
	 */
	@Override
	public void setId(String id) {
		this.id = id;
	}


    // **************************************************************
    // Attributable IMPLEMENTATION
    // **************************************************************

    /**
     * Sets an named attribute for this object.
     *
     * @param name  the name of the attribute.
     * @param value the new value of the attribute.
     */
    public void setAttribute(String name, Object value) {
        if (attributes == null)
            attributes = new Attributable.AttributeHelper();
        attributes.setAttribute(name, value);
    }

    /**
     * @param name the name of the attribute of interest.
     * @return an object representing the named attributed for this object.
     */
    public Object getAttribute(String name) {
        if (attributes == null)
            return null;
        else
            return attributes.getAttribute(name);
    }

    /**
     * @return an iterator of the attributes that this object has.
     */
    public Iterator<String> getAttributeNames() {
        if (attributes == null)
            return null;
        else
            return attributes.getAttributeNames();
    }

    
    // **************************************************************
    // static method
    // **************************************************************

	
	public static int[][] calculeteSPSArray(HaplotypeModel h1, HaplotypeModel h2){
		
		
		char[][] m1 = h1.getCharMatrix();
		char[][] m2 = h2.getCharMatrix();
		int[][] sps = calculateSPSCore(m1,m2);
		return sps;
		
	}
	
	public static int calculeteSPS(HaplotypeModel h1, HaplotypeModel h2){
		
		int sps = 0;
		char[][] m1 = h1.getCharMatrix();
		char[][] m2 = h2.getCharMatrix();
		int[][] spsArray = calculateSPSCore(m1, m2);
		for (int i = 0; i < spsArray.length; i++) {
			for (int j = 0; j < spsArray[i].length; j++) {
				sps += spsArray[i][j];
			}
		}
		return sps;
		
	}

	private static int[][] calculateSPSCore(char[][] m1, char[][] m2){
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

	

	public static Alignment swapAlignment(Alignment alignment){
		
		SimpleAlignment newAlignment = new SimpleAlignment();
		
		int seqCount = alignment.getSequenceCount();
		int siteCount = alignment.getSiteCount();
		
		for (int i = 0; i < seqCount; i++) {
			StringBuilder sb = new StringBuilder(siteCount);
			for (int j = 0; j < siteCount; j++) {
				int r = rand.nextInt(seqCount);
				char c = alignment.getAlignedSequenceString(r).charAt(j);
				sb.append(c);
			}
			
			Haplotype seq = new Haplotype(alignment.getTaxon(i), sb.toString());
			newAlignment.addSequence(seq);
		}
		return newAlignment;
	}

}