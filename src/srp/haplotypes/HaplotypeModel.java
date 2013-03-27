package srp.haplotypes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import dr.app.beauti.types.OperatorType;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel.TreeChangedEvent;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
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


public class HaplotypeModel extends AbstractHaplotypeModel  {

	public static final char GAP = '-';
	public static final String TAXON_PREFIX = "Hap_";
	
	private static final long serialVersionUID = -5057514703825711955L;

	private static Random rand = new Random();

//	int haplotypesCount;
	AlignmentMapping aMap;
	SwapInfo swapInfo = new SwapInfo();

	
	
	@Deprecated
	char[][] matrix;
	@Deprecated
	private int[] swapInfoOld = new int[NULL_SWAPINFO.length]; 
	//TODO: Expand to a class later, use this to speed up likelihood calculation
	@Deprecated
	public static final int[] NULL_SWAPINFO = new int[4];
	private static final String MODEL_NAME = "HaplotypeModel";
	
	
	
	
	
	
	
	
	
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
		haplotype.setDataType(getDataType());
	    haplotypes.add(haplotype);
	
	}

	public HaplotypeModel(AlignmentMapping aMap, int hapCount) {
		super(MODEL_NAME);
		setupData(aMap);
				
		initSeqs(hapCount);
		swapInfoOld = new int[4];//TODO remove later
		
	}

	public HaplotypeModel(AlignmentMapping aMap, Alignment trueAlignment) {
		super(MODEL_NAME);
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
		swapInfoOld[0] = hapIndex;
		swapInfoOld[1] = pos;
		swapInfoOld[2] = matrix[hapIndex][pos];
		swapInfoOld[3] = c;

		swapInfo.storeOperation(Operation.SWAPBASE, swapInfoOld);
		matrix[hapIndex][pos] = c;

		getHaplotype(hapIndex).setCharAt(pos, c);
		fireModelChanged();
		
	}
	public void fireModelChanged(){
//		for (TreeChangedEvent treeChangedEvent : treeChangedEvents) {
		listenerHelper.fireModelChanged(this);//, treeChangedEvent);
//		}
//		treeChangedEvents.clear();
	}


	@Deprecated
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
		swapInfoOld = swapInfo.getSwapInfo();
		matrix[swapInfoOld[0]][swapInfoOld[1]] = (char) swapInfoOld[2];
		getHaplotype(swapInfoOld[0]).setCharAt(swapInfoOld[1], (char)swapInfoOld[2]);
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
		return swapInfo.getSwapInfo();
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

	
	@Deprecated
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

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// TODO Auto-generated method stub
		System.err.println("Call handleModelChangedEvent");
		
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		// TODO Auto-generated method stub
		System.err.println("Call handleVariableChangedEvent");
	}

	@Override
	protected void storeState() {

//		System.err.println("Call storeState");
		swapInfo.storeOperation(Operation.NONE, null);
	}

	@Override
	protected void restoreState() {
//		System.err.println("Call restoreState");

		if(Operation.NONE != swapInfo.getOperation()){
			reject();
		}
	}

	@Override
	protected void acceptState() {
		//Do nothing
	}

}