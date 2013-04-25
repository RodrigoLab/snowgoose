package srp.haplotypes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.Iterator;

import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;

import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxon;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;
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
	public static final String TAXON_PREFIX = "hap_";
	
	private static final String MODEL_NAME = "HaplotypeModel";
	private static final long serialVersionUID = -5057514703825711955L;



//	int haplotypesCount;
	AlignmentMapping aMap;
	private SwapInfo swapInfo = new SwapInfo();
	private boolean isEdit;

	
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
			randomSeq(i);

		}
	}
	
	private void addHaplotype(Haplotype haplotype) {
		haplotype.setDataType(getDataType());
	    haplotypes.add(haplotype);
	
	}

	public HaplotypeModel(AlignmentMapping aMap, int hapCount) {
		super(MODEL_NAME);
		setupData(aMap);//		swapBase(hapIndex, pos);
		
		initSeqs(hapCount);
	}

	public HaplotypeModel(AlignmentMapping aMap, Alignment trueAlignment) {
		super(MODEL_NAME);
		setupData(aMap);
		setupAlignment(trueAlignment);
	
	}

	
	public void randomSeq(int hapIndex) {
		
		for (int i = 0; i < getHaplotypeLength(); i++) {
//			char newChar = GAP;
//			int size = aMap.mapToSrp[i].size();
//			if (size != 0) {
//				int srpIndex = aMap.mapToSrp[i].get(rand.nextInt(0, size));
//				newChar = aMap.getShortReadCharAt(srpIndex, i);
//			}
			char newChar = (char) aMap.nextBaseAt(i);
//			int[] swapInfoArray = swapHaplotypeBase(hapIndex, posChar);

			getHaplotype(hapIndex).setCharAt(i, newChar);

		}
	}
	public void startHaplotypeOperation(){
		isEdit = true;
	}
	public void endHaplotypeOperation(){
		isEdit = false;
		fireModelChanged();
	}
	
	public int[] getNextBase(){
		return aMap.nextBase();
	}
	public int[] getNextBaseUniform(){
		return aMap.nextBaseUniform();
	}
	
	public void storeOperationRecord(Operation op, Object opRecord){
		swapInfo.storeOperation(op, opRecord);
	}
	public int[] swapHaplotypeBase(int hapIndex, int[] posChar){
//		replaceHaplotypeCharAt(hapIndex, posChar);
		
		int[] swapRecord = new int[4];
		swapRecord[0] = hapIndex;
		swapRecord[1] = posChar[0];
		swapRecord[2] = posChar[1];
		swapRecord[3] = swapHaplotypeCharAt(hapIndex, posChar[0], posChar[1]);

//		System.out.println(Arrays.toString(swapRecord));
//		System.out.println(getHaplotype(hapIndex).charAt(posChar[0]));

//		swapInfoArray[2] = replaceHaplotypeCharAt(hapIndex, posChar);
		return swapRecord;
	}

	private int swapHaplotypeCharAt(int hapIndex, int pos, int newChar){
		int oldChar = getHaplotype(hapIndex).getChar(pos); 
		getHaplotype(hapIndex).setCharAt(pos, (char) newChar);
		return oldChar;
		
	}
	private void resetHaplotypeToOldChar(int[] swapArray){
		swapHaplotypeCharAt(swapArray[0], swapArray[1], swapArray[3]);
	}
	
//	private int swapHaplotypeCharAt(int hapIndex, int[] posChar){
//		int oldChar = getHaplotype(hapIndex).getChar(posChar[0]); 
//		getHaplotype(hapIndex).setCharAt(posChar[0], (char) posChar[1]);
//		return oldChar;
//		
//	}
//	
//	

	public void reject() {//FIXME for different swapInfo/Operation
		Operation op = swapInfo.getOperation();
		int[] temp;
		switch (op) {
			case SWAPBASE: 
				
				temp = swapInfo.getSwapInfoSWAPBASE();
				resetHaplotypeToOldChar(temp);
				break;
			case UNIFORMSWAPBASE:
				
				temp = swapInfo.getSwapInfoSWAPBASE();
				resetHaplotypeToOldChar(temp);
				break;
			case SWAPMULTI:
				Deque<int[]> swapMulti = swapInfo.getSwapInfoSWAPMULTI();
				
				for (Iterator<int[]> iterator = swapMulti.descendingIterator(); iterator
						.hasNext();) {

					temp= iterator.next();
					resetHaplotypeToOldChar(temp);
				}

				break;
			case SWAPSECTION:
				temp = swapInfo.getSwapHaplotypeRecord();
				for (int i = 0; i < temp.length; i++) {
					haplotypes.get(temp[i]).restoreState();
				}
				break;
			case RECOMB:
				temp = swapInfo.getSwapHaplotypeRecord();
				for (int i = 0; i < temp.length; i++) {
					haplotypes.get(temp[i]).restoreState();
				}
				break;
				
				
			default:
				throw new IllegalArgumentException("Unknown operation type: "+op);
				
		}
		
	}
	


	@Override
	public void fireModelChanged(){
//		for (TreeChangedEvent treeChangedEvent : treeChangedEvents) {
		listenerHelper.fireModelChanged(this);//, treeChangedEvent);
//		}
//		treeChangedEvents.clear();
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
	
	public Operation getOperation() {
		 
		return swapInfo.getOperation();
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



    
	
	public String[] toStringArray() {
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
//		System.err.println("Call restoreState - haplotypeModel");

		if(Operation.NONE != swapInfo.getOperation()){
			reject();
		}
	}

	@Override
	protected void acceptState() {
		//Do nothing
	}

	public AlignmentMapping getAlignmentMapping() {
		return aMap;
	}


}