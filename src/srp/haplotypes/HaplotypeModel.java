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
	
	private static final String MODEL_NAME = "HaplotypeModel";
	private static final long serialVersionUID = -5057514703825711955L;

	private static Random rand = new Random();

//	int haplotypesCount;
	AlignmentMapping aMap;
	SwapInfo swapInfo = new SwapInfo();

	
	
	
	
	
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
		
    }
	
	private void initSeqs(int hapCount){
		
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

		
	}

	public HaplotypeModel(AlignmentMapping aMap, Alignment trueAlignment) {
		super(MODEL_NAME);
		setupData(aMap);
		setupAlignment(trueAlignment);

	
	}

	
	public void randomSeq(int hapIndex) {
		randomSeq(hapIndex, 0, getHaplotypeLength());
		
	}
	
	
	public void randomSeq(int hapIndex, int start, int end){
		
		for (int p = start; p < end; p++) {
			swapBase(hapIndex, p);

		}
		
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

	private void swapBase(int hapIndex, int pos, char newChar){
		int[] swapInfo2 = new int[4];
		swapInfo2[0] = hapIndex;
		swapInfo2[1] = pos;
		swapInfo2[2] = getHaplotype(hapIndex).getChar(pos); //matrix[hapIndex][pos];
		swapInfo2[3] = newChar;

		swapInfo.storeOperation(Operation.SWAPBASE, swapInfo2);
		
		getHaplotype(hapIndex).setCharAt(pos, newChar);
	
		
		fireModelChanged();
		
	}
	public void fireModelChanged(){
//		for (TreeChangedEvent treeChangedEvent : treeChangedEvents) {
		listenerHelper.fireModelChanged(this);//, treeChangedEvent);
//		}
//		treeChangedEvents.clear();
	}



	public void reject() {//FIXME for different swapInfo/Operation
		int[] swapInfo2 = swapInfo.getSwapInfoIntArray();
		
		getHaplotype(swapInfo2[0]).setCharAt(swapInfo2[1], (char)swapInfo2[2]);
	}


	
	public int calculateSPS(){
		int sps = 0;
		
		for (int i = 0; i < getHaplotypeCount(); i++) {
			Haplotype h1 = getHaplotype(i);
			for (int j = 0; j < i; j++) {
				Haplotype h2 = getHaplotype(j);
				for (int pos = 0; pos < getHaplotypeLength(); pos++) {
					int c = h1.charAt(pos) - h2.charAt(pos);
					sps += (c==0)? 0: 1;
				}
			}
		}
		return sps;
	}

	public SwapInfo getSwapInfo() {
		return swapInfo;
	}

	
	@Deprecated
	public int[] getSwapInfoIntArray() {
		return swapInfo.getSwapInfoIntArray();
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
		
		int[][] sps = calculateSPSCore(h1,h2);
		return sps;
		
	}
	
	public static int calculeteSPS(HaplotypeModel h1, HaplotypeModel h2){
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
	
	private static int[][] calculateSPSCore(HaplotypeModel h1, HaplotypeModel h2){
			int hapCount = h1.getHaplotypeCount();
			int seqLength = h1.getHaplotypeLength();
			int sps[][] = new int[hapCount][hapCount];
			if (seqLength != h2.getHaplotypeLength()){
				System.err.println("Incompariable alignments lenght: "+seqLength +" and "+  h2.getHaplotypeLength());
			}
			if (hapCount != h2.getHaplotypeCount()){
				System.err.println("Different number of haplotypes: "+hapCount +" and "+  h2.getHaplotypeCount());
			}
			String[] s1 = h1.toStringArray();
			String[] s2 = h2.toStringArray();
			for (int i = 0; i < hapCount; i++) {
				for (int j = 0; j < hapCount; j++) {
					sps[i][j] = caluclateSPSSingle(s1[i], s2[j]);
				}
			}
			return sps;
			
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
	
	private String[] toStringArray() {
		String[] string = new String[getHaplotypeCount()];
		for (int i = 0; i < string.length; i++) {
			string[i] = getHaplotypeString(i);
		}
		return string;
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