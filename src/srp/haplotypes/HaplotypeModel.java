package srp.haplotypes;


import java.util.ArrayList;
import java.util.Arrays;

import srp.haplotypes.SwapInfo.Operation;

import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxon;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
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

	private static final char[] VALID_CHARS = new char[4];
	
	public static final char GAP = '-';
	public static final String TAXON_PREFIX = "hap_";
	
	private static final String MODEL_NAME = "HaplotypeModel";
	private static final long serialVersionUID = -5057514703825711955L;

	private static final int INDEX_OF_LAST_VALID_CHARS = 3;

	

//	int haplotypesCount;
	AlignmentMapping aMap;
	private SwapInfo swapInfo = new SwapInfo();
	
	
	private boolean isEdit;
	
	
	private HaplotypeModel(AlignmentMapping aMap){
		super(MODEL_NAME);

		this.aMap = aMap;
		haplotypesLength = this.aMap.getLength();
		
		haplotypes = new ArrayList<Haplotype>();
		setDataType(Nucleotides.INSTANCE);

		System.arraycopy(Nucleotides.NUCLEOTIDE_CHARS, 0, VALID_CHARS, 0, VALID_CHARS.length);
		storedCumSumFrequency[INDEX_OF_LAST_VALID_CHARS]=1;
		
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

//
//
//		int i=0;
//		Taxon t = new Taxon(TAXON_PREFIX+i); 
//		Haplotype haplotype = new Haplotype(t, tempSeq);
//		addHaplotype(haplotype);
//		randomHaplotype(i);
//
//		for (i = 1; i < hapCount; i++) {
//			t = new Taxon(TAXON_PREFIX+i); 
//			Haplotype haplotype2 = new Haplotype(t, haplotype.getSequenceString());
//			addHaplotype(haplotype2);
//		}


	}
	

	private void randomHaplotype(int hapIndex) {
		Haplotype haplotype = haplotypes.get(hapIndex);
		for (int i = 0; i < haplotypesLength; i++) {
			//TODO: Which is the "good" starting point
//			char newChar = (char) aMap.nextBaseAt(i);
			char newChar = (char) aMap.nextBaseEqualFreq();
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
	

	public void swapHaplotypeColumn(int[] posChar){

		int[] allOldChars = new int[getHaplotypeCount()];
			
		for (int i = 0; i < getHaplotypeCount(); i++) {
			Haplotype haplotype = haplotypes.get(i);
			allOldChars[i] = haplotype.replaceCharAt(posChar[0], posChar[1]);
		}
			
		storeOperationRecord(Operation.SWAPCOLUMN, posChar, allOldChars);


	}
	public int[] swapHaplotypeSingleBase(Operation op, int[] posChar){
		
		int[] swapInfoArray = new int[4];
		swapInfoArray[0] = MathUtils.nextInt(getHaplotypeCount());
		swapInfoArray[1] = posChar[0];
		swapInfoArray[2] = posChar[1];
		
		Haplotype haplotype = haplotypes.get(swapInfoArray[0]);
		swapInfoArray[3] = haplotype.replaceCharAt(posChar[0], posChar[1]);
		
		storeOperationRecord(op, swapInfoArray);
		return swapInfoArray;
//		swapInfo.storeOperation(op, swapInfoArray);

	}

	@Deprecated
	public void swapHaplotypeMultiBases(Operation op, int[] allNewChars){
		int hapIndex = MathUtils.nextInt(getHaplotypeCount());
		
		Haplotype haplotype = haplotypes.get(hapIndex);
		
		int[] allOldChars = new int[allNewChars.length];
		for (int i = 0; i < haplotypesLength; i++) {
			int newChar = allNewChars[i];
			if(newChar>0){
				allOldChars[i] = haplotype.replaceCharAt(i, newChar);
			}
		}
		storeOperationRecord(op, new int[]{hapIndex}, allNewChars, allOldChars);
//		swapInfo.storeOperation(op, new int[]{hapIndex}, allPosChars);
		
	}
	public void swapHaplotypeMultiBases(Operation op, int[][] allPosChars){
		int hapIndex = MathUtils.nextInt(getHaplotypeCount());
		
		Haplotype haplotype = haplotypes.get(hapIndex);

		for (int i = 0; i < haplotypesLength; i++) {
			int newChar = allPosChars[0][i];
			if(newChar>0){
				allPosChars[1][i] = haplotype.replaceCharAt(i, newChar);
			}
		}
		storeOperationRecord(op, hapIndex, allPosChars);
//		storeOperationRecord(op, new int[] {hapIndex}, allPosChars[0], allPosChars[1]);
//		swapInfo.storeOperation(op, new int[]{hapIndex}, allPosChars);
		
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
	private void resetHaplotypeToOldChar(int[] swapRecord){
		Haplotype haplotype = haplotypes.get(swapRecord[0]);
		haplotype.replaceCharAt(swapRecord[1], swapRecord[3]);

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
//	public Operation getOperation() {
//		return swapInfo.getOperation();
//	}
	public void storeOperationRecord(Operation op, int[]... opRecord){
		swapInfo.storeOperation(op, opRecord);
	}
	
	private void storeOperationRecord(Operation op, int hapIndex,
			int[][] opRecord) {
		swapInfo.storeOperation(op, hapIndex, opRecord);		
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
//		int[] temp;
		switch (op) {
		
		case NONE:
			break;
		case SWAPSINGLE:

			int[] temp = swapInfo.getSwapInfoSWAPBASE();
			//			long time1 = System.currentTimeMillis();
//			for (int i = 0; i < 1e9; i++) {
			resetHaplotypeToOldChar(temp);
//		replaceHaplotypeCharAt(haplotypes.get(temp[0]), temp[2], temp[3]);
//			}
//			long time2 = System.currentTimeMillis();
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
			int[] swapSection = swapInfo.getSwapInfoSWAPSECTION();
			
//			time1 = System.currentTimeMillis();
////			for (int i = 0; i < 1e8; i++) {
//				haplotypes.get(temp[0]).restoreState();
////			}
//			time2 = System.currentTimeMillis();

//			System.out.println("Multi: "+(time2 - time1) + "\t");

			
			
			for (int i = 0; i < 2; i++) {
				haplotypes.get(swapSection[i]).restoreState();
			}
			break;
		case SWAPCOLUMN:
			int[][] swapColumn = swapInfo.getSwapInfoSWAPCOLUMN();
			int pos = swapColumn[0][0];
			int[] allOldChars2 = swapColumn[1];

			for (int i = 0; i < getHaplotypeCount(); i++) {
				haplotype = haplotypes.get(i);
				haplotype.replaceCharAt(pos, allOldChars2[i]);
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
	
	public double getLogqFrequency(int oldChar, int newChar){
		return storedLogqMatrix[oldChar][newChar];
	}
	
	public int[] getNextBaseFrequency(Parameter frequency) {
	
		checkFrequencyParameter(frequency);

		int[] tempPosChar = new int[2];
		tempPosChar[0] = MathUtils.nextInt(getHaplotypeLength());

		double d = MathUtils.nextDouble();

		for (int i = 0; i < INDEX_OF_LAST_VALID_CHARS; i++) {
			if (d <= storedCumSumFrequency[i]) {
				tempPosChar[1] = VALID_CHARS[i];
				return tempPosChar;
			}
		}
		tempPosChar[1] = VALID_CHARS[INDEX_OF_LAST_VALID_CHARS];
		return tempPosChar;
	}
	

	private void checkFrequencyParameter(Parameter frequency) {

		for (int i = 0; i < storedFrequency.length; i++) {
			if(storedFrequency[i]!= frequency.getParameterValue(i)){

				for (int j = i; j < storedFrequency.length; j++) {
					storedFrequency[j] = frequency.getParameterValue(j);
					logFreq[j] = Math.log(storedFrequency[j]);
				}			
				storedCumSumFrequency[0] = storedFrequency[0];
				storedCumSumFrequency[1] = storedCumSumFrequency[0]+storedFrequency[1];
				storedCumSumFrequency[2] = storedCumSumFrequency[1]+storedFrequency[2];

				for (int j = 0; j < logFreq.length; j++) {
					for (int k = j+1; k < logFreq.length; k++) {
						storedLogqMatrix[j][k] = logFreq[j]-logFreq[k];
						storedLogqMatrix[k][j] = -storedLogqMatrix[j][k];
					}
//					System.out.println(Arrays.toString(storedLogqMatrix[j]));
				}
				
				break;
			}
		}
	}

	private void checkFrequencyParameterTune(Parameter frequency) {

		for (int i = 0; i < storedFrequency.length; i++) {
			if(storedFrequency[i]!= frequency.getParameterValue(i)){
				System.out.println(i);
long time1 = System.currentTimeMillis();
for (int t = 0; t < 1e7; t++) {
	
				for (int j = i; j < storedCumSumFrequency.length; j++) {
					storedFrequency[j] = frequency.getParameterValue(j);
					logFreq[j] = Math.log(storedFrequency[j]);
				}			
//				storedFrequency = frequency.getParameterValues();
				storedCumSumFrequency[0] = storedFrequency[0];
				storedCumSumFrequency[1] = storedCumSumFrequency[0]+storedFrequency[1];
				storedCumSumFrequency[2] = storedCumSumFrequency[1]+storedFrequency[2];
				// Too short for a loop?
				// for (int j = 1; j < INDEX_OF_LAST_VALID_CHARS; j++) {
				// storedCumSumFrequency[j] = storedCumSumFrequency[j - 1]
				// + storedFrequency[j];
				// }
				
//				double[] logFreq = new double[4];
//				for (int j = 0; j < storedCumSumFrequency.length; j++) {
//					logFreq[j] = Math.log(storedFrequency[j]);
//				}
//				System.err.println(Arrays.toString(storedFrequency));
				for (int j = 0; j < storedFrequency.length; j++) {
					for (int k = j+1; k < storedFrequency.length; k++) {
//						storedLogqMatrix[j][k] = Math.log(storedFrequency[j]/storedFrequency[k]);
						storedLogqMatrix[j][k] = logFreq[j]-logFreq[k];
						storedLogqMatrix[k][j] = -storedLogqMatrix[j][k];
					}
//					System.out.println(Arrays.toString(storedLogqMatrix[j]));
				}
//				System.out.println();
//				for (int j = 0; j < storedFrequency.length; j++) {
//					for (int k = 0; k < storedFrequency.length; k++) {
////						storedLogqMatrix[j][k] = Math.log(storedFrequency[j]/storedFrequency[k]);
//						storedLogqMatrix[j][k] = logFreq[j]-logFreq[k];
//					}
////					System.out.println(Arrays.toString(storedLogqMatrix[j]));
//				}

//				System.out.println(i +"\t"+ Arrays.toString(storedFrequency));
//				System.err.println(Arrays.toString(storedCumSumFrequency));
}
long time2 = System.currentTimeMillis();

System.out.println((time2 - time1) + "\t");
				
				break;
			}
		}
	}
	
	private double[] logFreq = new double[4];
	private double[] storedFrequency = new double[4];
	private double[] storedCumSumFrequency = new double[4];
	private double[][] storedLogqMatrix = new double[4][4];
}