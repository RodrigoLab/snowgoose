package srp.haplotypes;

import java.util.ArrayList;
import java.util.Arrays;

import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxon;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
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
	

	
	private HaplotypeModel(AlignmentMapping aMap){
			super(MODEL_NAME);
	//		setupData(aMap);
			this.aMap = aMap;
			this.haplotypesLength = this.aMap.getLength();
			
			haplotypes = new ArrayList<Haplotype>();
			setDataType(Nucleotides.INSTANCE);
			
		}

	private void addHaplotype(Haplotype haplotype) {
		haplotype.setDataType(getDataType());
	    haplotypes.add(haplotype);
	
	}

	private void initSeqs(int hapCount){
		
		char[] temp = new char[haplotypesLength];
		Arrays.fill(temp, DataType.GAP_CHARACTER);
		String tempSeq = String.valueOf(temp);
		
		for (int i = 0; i < hapCount; i++) {
			Taxon t = new Taxon(TAXON_PREFIX+i); 
			Haplotype haplotype = new Haplotype(t, tempSeq);
			addHaplotype(haplotype);
			randomHaplotype(i);

		}

	}
	

	private void randomHaplotype(int hapIndex) {
		Haplotype haplotype = haplotypes.get(hapIndex);
		for (int i = 0; i < haplotypesLength; i++) {
			char newChar = (char) aMap.nextBaseAt(i);
			haplotype.setCharAt(i, newChar);
	
		}
	}
	
	public HaplotypeModel(AlignmentMapping aMap, int hapCount) {
		this(aMap);		
		initSeqs(hapCount);
		
	}

	public HaplotypeModel(AlignmentMapping aMap, Alignment trueAlignment) {
		this(aMap);

		for (int i = 0; i < trueAlignment.getSequenceCount(); i++) {
			Haplotype haplotype = new Haplotype(trueAlignment.getSequence(i));
			addHaplotype(haplotype);
		}

	}	
	

	
	public int[] swapHaplotypeSingleBase(int hapIndex, int[] posChar){
		
		int[] swapRecord = new int[4];
		swapRecord[0] = hapIndex;
		swapRecord[1] = posChar[0];
		swapRecord[2] = posChar[1];
//		swapRecord[3] = swapHaplotypeCharAt(hapIndex, posChar[0], posChar[1]);
		
		Haplotype haplotype = haplotypes.get(hapIndex);
		swapRecord[3] = haplotype.replaceCharAt(posChar[0], posChar[1]);

//		swapRecord[3] = haplotype.getChar(posChar[0]);
//		setCharAt(pos, (char) newChar);
		
		
//		System.out.println(Arrays.toString(swapRecord));
//		System.out.println(getHaplotype(hapIndex).charAt(posChar[0]));

//		swapInfoArray[2] = replaceHaplotypeCharAt(hapIndex, posChar);
		return swapRecord;
	}


	public int[] swapHaplotypeMultiBases(int hapIndex, int[] allPosChars, int[] oldChars){
		
		Haplotype haplotype = haplotypes.get(hapIndex);

		for (int i = 0; i < haplotypesLength; i++) {
			int newChar = allPosChars[i];
			if(newChar>0){
				oldChars[i] = haplotype.replaceCharAt(i, newChar);
			}
		}
		return oldChars;
	}

	
	
//	private int replaceHaplotypeCharAt(int hapIndex, int pos, int newChar){
//		
//		int oldChar = haplotypes.get(hapIndex).replaceCharAt(pos, (char) newChar);
//		return oldChar;
//		
//	}
//	private static int replaceHaplotypeCharAt(Haplotype haplotype, int pos, int newChar){
//		int oldChar = haplotype.replaceCharAt(pos, (char) newChar);
//		return oldChar;
//		
//	}
	private void resetHaplotypeToOldChar(int[] swapArray){
		Haplotype haplotype = haplotypes.get(swapArray[0]);
		haplotype.replaceCharAt(swapArray[1], swapArray[3]);

//		swapHaplotypeCharAt(swapArray[0], swapArray[1], swapArray[3]);
	}
	
//	private int swapHaplotypeCharAt(int hapIndex, int[] posChar){
//		int oldChar = getHaplotype(hapIndex).getChar(posChar[0]); 
//		getHaplotype(hapIndex).setCharAt(posChar[0], (char) posChar[1]);
//		return oldChar;
//		
//	}
//	
//	

	public SwapInfo getSwapInfo() {
		return swapInfo;
	}
	public Operation getOperation() {
		return swapInfo.getOperation();
	}
	public void storeOperationRecord(Operation op, Object... opRecord){
		swapInfo.storeOperation(op, opRecord);
	}



	
	public int calculateSPS(){
		int sps = 0;
		
		for (int i = 0; i < getHaplotypeCount(); i++) {
			Haplotype h1 = haplotypes.get(i);
			for (int j = 0; j < i; j++) {
				Haplotype h2 = haplotypes.get(j);
				for (int pos = 0; pos < haplotypesLength; pos++) {
					int c = h1.getChar(pos) - h2.getChar(pos);
					sps += (c==0)? 0: 1;
				}
			}
		}
		return sps;
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



    
	
	@Override
		public void fireModelChanged(){
	//		for (TreeChangedEvent treeChangedEvent : treeChangedEvents) {
			listenerHelper.fireModelChanged(this);//, treeChangedEvent);
	//		}
	//		treeChangedEvents.clear();
		}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		System.err.println("Call handleModelChangedEvent");
		
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		System.err.println("Call handleVariableChangedEvent");
	}

	@Override
	protected void storeState() {

		// System.err.println("Call storeState");
		swapInfo.storeOperation(Operation.NONE);
	}

	@Override
	protected void restoreState() {
		// System.err.println("Call restoreState - haplotypeModel");
		reject();
		swapInfo.storeOperation(Operation.NONE);
	}

	//	private int swapHaplotypeCharAt(int hapIndex, int[] posChar){
	//		int oldChar = getHaplotype(hapIndex).getChar(posChar[0]); 
	//		getHaplotype(hapIndex).setCharAt(posChar[0], (char) posChar[1]);
	//		return oldChar;
	//		
	//	}
	//	
	//	

	public void reject() {
		Operation op = swapInfo.getOperation();
		int[] temp;
		switch (op) {
		case NONE:
			break;
		case SWAPSINGLE:

			temp = swapInfo.getSwapInfoSWAPBASE();
			long time1 = System.currentTimeMillis();
//			for (int i = 0; i < 1e9; i++) {
		resetHaplotypeToOldChar(temp);
//		replaceHaplotypeCharAt(haplotypes.get(temp[0]), temp[2], temp[3]);
//			}
			long time2 = System.currentTimeMillis();
//
//			System.out.println("Single: "+(time2 - time1) + "\t");
//			System.exit(0);
			break;

		case SWAPMULTI:
			
			int hapIndex = swapInfo.getHapIndex();
			int[][] swapMulti = swapInfo.getSwapInfoSWAPMULTI();
//			int[] allNewChars = swapMulti[0];
			int[] allOldChars = swapMulti[1];
			
			
			Haplotype haplotype = haplotypes.get(hapIndex);
			for (int i = 0; i < allOldChars.length; i++) {
				int oldChar = allOldChars[i];
				if(oldChar>0){
					haplotype.replaceCharAt(i, oldChar);
				}
			}
			break;
		
		case SWAPSECTION:
			temp = swapInfo.getSwapHaplotypeRecord();
			
//			time1 = System.currentTimeMillis();
////			for (int i = 0; i < 1e8; i++) {
//				haplotypes.get(temp[0]).restoreState();
////			}
//			time2 = System.currentTimeMillis();

//			System.out.println("Multi: "+(time2 - time1) + "\t");

			
			
			for (int i = 0; i < 2; i++) {
				haplotypes.get(temp[i]).restoreState();
			}
			break;

		default:
			throw new IllegalArgumentException("Unknown operation type: " + op);

		}

	}

	@Override
	protected void acceptState() {
		//Do nothing
	}

	public void startHaplotypeOperation(){
		isEdit = true;
	}

	public void endHaplotypeOperation(){
		isEdit = false;
		fireModelChanged();
	}


	
	public AlignmentMapping getAlignmentMapping() {
		return aMap;
	}

	public char getHaplotypeCharAt(int hapIndex, int charIndex) {
		return haplotypes.get(hapIndex).getChar(charIndex);
	}
	
	

	@Deprecated
	public String[] toStringArray() {
		String[] string = new String[getHaplotypeCount()];
		for (int i = 0; i < string.length; i++) {
			string[i] = getHaplotypeString(i);
		}
		return string;
	}

}