package dr.rj;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jebl.evolution.sequences.Sequence;
import srp.evolution.haplotypes.old.OldAbstractHaplotypeModel;
import srp.evolution.haplotypes.old.OldHapOperation;
import srp.evolution.haplotypes.old.OldHapSwapInfo;
import srp.evolution.haplotypes.old.OldHaplotype;
import srp.evolution.shortreads.AlignmentMapping;
import srp.haplotypes.AbstractHaplotypeModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.ext.SeqGenExt;
import dr.ext.TreeLikelihoodExt;
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


public class SSHaplotypeModel extends OldAbstractHaplotypeModel  {
	
	private static final long serialVersionUID = 7277181804036354716L;
	private static final String MODEL_NAME = "HaplotypeModel";

	
	private static final int HAP_INDEX = OldHapSwapInfo.SWAPBASE_HAP_INDEX;
	private static final int POS_INDEX = OldHapSwapInfo.SWAPBASE_POS_INDEX;
	private static final int NEW_CHAR_INDEX = OldHapSwapInfo.SWAPBASE_NEW_CHAR_INDEX;
	private static final int OLD_CHAR_INDEX = OldHapSwapInfo.SWAPBASE_OLD_CHAR_INDEX;

	private static final int NUCLEOTIDE_STATES[] = Nucleotides.NUCLEOTIDE_STATES;
	private static final char[] VALID_CHARS = initValidChars4();
	private static final char[] initValidChars4() {
		char[] validChar = new char[4];
		System.arraycopy(Nucleotides.NUCLEOTIDE_CHARS, 0, validChar, 0, validChar.length);
		return validChar;
	}
	private static final int INDEX_OF_LAST_VALID_CHARS = VALID_CHARS.length-1;
	
//	private static final DataType DATA_TYPE = Nucleotides.INSTANCE;
	public static final char GAP = '-';
	public static final String TAXON_PREFIX = "hap_";
	
	AlignmentMapping aMap;
	private OldHapSwapInfo swapInfo = new OldHapSwapInfo();
	
	private boolean isEdit;
	int[] swapBaseRecord = new int[4];
	
	//
	// SSModel
	//
	private static final int SCALE_TAXA_FACTOR = 2;
	private int trueHaplotypeCount;
	private int initHaplotypeCount = 16;
	


	private SSHaplotypeModel(AlignmentMapping aMap){
		super(MODEL_NAME);

		this.aMap = aMap;
		haplotypesLength = this.aMap.getLength();
		
		haplotypes = new ArrayList<OldHaplotype>();
		setDataType(Nucleotides.INSTANCE);
		
//		storedCumSumFrequency[INDEX_OF_LAST_VALID_CHARS]=1;
	}



	public SSHaplotypeModel(Alignment shortReads, int hapCount) {
		this(new AlignmentMapping(shortReads), hapCount);
	}
	
	public SSHaplotypeModel(AlignmentMapping aMap, int hapCount) {
		this(aMap);		
		initSeqs(hapCount);
		appendEmptySequences();
	}
	
	public SSHaplotypeModel(AlignmentMapping aMap, Alignment trueAlignment) {
		this(aMap);
		for (int i = 0; i < trueAlignment.getSequenceCount(); i++) {
			OldHaplotype haplotype = new OldHaplotype(trueAlignment.getSequence(i));
			addHaplotype(haplotype);
		}
		appendEmptySequences();
	}
	
	private void initSeqs(int hapCount){
		
		char[] temp = new char[haplotypesLength];
		Arrays.fill(temp, DataType.GAP_CHARACTER);
		String tempSeq = String.valueOf(temp);
		
		for (int i = 0; i < hapCount; i++) {
			Taxon t = new Taxon(TAXON_PREFIX+i); 
			OldHaplotype haplotype = new OldHaplotype(t, tempSeq);
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
		OldHaplotype haplotype = haplotypes.get(hapIndex);
		for (int i = 0; i < haplotypesLength; i++) {
			//TODO: Which is the "good" starting point
//			char newChar = (char) aMap.nextBaseAt(i);
			char newChar = (char) aMap.nextBaseEqualFreq();
			haplotype.setCharAt(i, newChar);
			
		}
	}
	private void addHaplotype(OldHaplotype haplotype) {
		haplotype.setDataType(getDataType());
	    haplotypes.add(haplotype);
	
	}
	private void appendEmptySequences() {
		
		trueHaplotypeCount = getHaplotypeCount();
		int haplotypeCount = getHaplotypeCount();

		haplotypeCount = initHaplotypeCount;
		while(trueHaplotypeCount > haplotypeCount ){
			haplotypeCount *= SCALE_TAXA_FACTOR;
		}
		
		char[] temp = new char[haplotypesLength];
		Arrays.fill(temp, DataType.GAP_CHARACTER);
		String tempSeq = String.valueOf(temp);
		
		for (int i = trueHaplotypeCount; i < haplotypeCount; i++) {
			Taxon t = new Taxon(TAXON_PREFIX+i); 
			OldHaplotype haplotype = new OldHaplotype(t, tempSeq);
			addHaplotype(haplotype);
//			randomHaplotype(i);
		}

		
	}
	
	private char setNewCharFromFrequency(){
		
		double d = MathUtils.nextDouble();
		
		for (int i = 0; i < INDEX_OF_LAST_VALID_CHARS; i++) {
			if (d <= storedCumSumFrequency[i]) {
				return VALID_CHARS[i];
			}
		}
		return VALID_CHARS[INDEX_OF_LAST_VALID_CHARS];
	}
	
	public int getTrueHaplotypeCount() {
		
		return trueHaplotypeCount;
	}



	//	private int swapHaplotypeCharAt(int hapIndex, int[] posChar){
	//		int oldChar = getHaplotype(hapIndex).getChar(posChar[0]); 
	//		getHaplotype(hapIndex).setCharAt(posChar[0], (char) posChar[1]);
	//		return oldChar;
	//		
	//	}
	//	
	//	

	public OldHapSwapInfo getSwapInfo() {
		return swapInfo;
	}



	public int[] getNextBaseFrequency(Parameter frequency) {
	
		checkFrequencyParameter(frequency);
	
		int[] tempPosChar = new int[2];
		tempPosChar[0] = MathUtils.nextInt(getHaplotypeLength());
	
		double d = MathUtils.nextDouble();
	
//		for (int i = 0; i < INDEX_OF_LAST_VALID_CHARS; i++) {
//			if (d <= storedCumSumFrequency[i]) {
//				tempPosChar[1] = VALID_CHARS[i];
//				return tempPosChar;
//			}
//		}
		tempPosChar[1] = setNewCharFromFrequency();
		return tempPosChar;
	}

	
	
//	private int swapHaplotypeCharAt(int hapIndex, int[] posChar){
//		int oldChar = getHaplotype(hapIndex).getChar(posChar[0]); 
//		getHaplotype(hapIndex).setCharAt(posChar[0], (char) posChar[1]);
//		return oldChar;
//		
//	}
//	
//	

	public int calculateSPS(){
		int sps = 0;
		
		for (int i = 0; i < getHaplotypeCount(); i++) {
			OldHaplotype h1 = haplotypes.get(i);
			for (int j = 0; j < i; j++) {
				OldHaplotype h2 = haplotypes.get(j);
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
		swapInfo.storeOperation(OldHapOperation.NONE);
	}

	@Override
	protected void restoreState() {
		// System.err.println("Call restoreState - haplotypeModel");
		reject();
		swapInfo.storeOperation(OldHapOperation.NONE);
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
		
		OldHapOperation op = swapInfo.getOperation();
//		int[] temp;
		switch (op) {
		
		case NONE:
			break;
		case PASS:
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
			
			
			OldHaplotype haplotype = haplotypes.get(hapIndex);
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
////////////////////////////////////////////////////////////////////////
	public double getLogqFrequency(int oldChar, int newChar){
		return storedLogqMatrix[NUCLEOTIDE_STATES[oldChar]][NUCLEOTIDE_STATES[newChar]];
	}
	
	public double getLogqFrequencyStates(int oldState, int newState){
		return storedLogqMatrix[oldState][newState];
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
//				storedCumSumFrequency[2] = storedCumSumFrequency[1]+storedFrequency[2];

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
	private double[] storedCumSumFrequency = new double[INDEX_OF_LAST_VALID_CHARS];
	private double[][] storedLogqMatrix = new double[4][4];


	public void simulateSequence(TreeLikelihoodExt treeLikelihood) {

        double substitutionRate = 1000.1;//(getHaplotypeCount()*getHaplotypeLength()) ;
        double damageRate = 0;
        SiteModel siteModel = treeLikelihood.getSiteModel();
        SubstitutionModel substitutionModel = siteModel.getSubstitutionModel();
        
        int[] initialSequence = aMap.getConsensusSequenceState();
        StringBuffer buffer = new StringBuffer();
        for (int i = 0; i < initialSequence.length; i++) {
            buffer.append(Nucleotides.INSTANCE.getChar(    initialSequence[i] ));
        }
        treeLikelihood.makeDirty();
        treeLikelihood.getLogLikelihood();
        System.err.println(buffer.toString());
//        sequence.setDataType(dataType);


//    	Arrays.fill(initialSequence, 1);
        SeqGenExt seqGen = new SeqGenExt(initialSequence, 
                substitutionRate, substitutionModel, siteModel,
                damageRate);
        
        Tree tree = treeLikelihood.getTreeModel();
        jebl.evolution.alignments.Alignment jeblAlignment = seqGen.simulate(tree);
        List<Sequence> sequenceList = jeblAlignment.getSequenceList();
        for (int j = 0; j < sequenceList.size(); j++) {
			

//        	System.err.println(getHaplotypeString(j));
			System.out.println(sequenceList.get(j).getString());
			OldHaplotype haplotype = getHaplotype(j);
			haplotype.setSequenceString(sequenceList.get(j).getString());
			System.out.println(getHaplotypeString(j));
			System.out.println();
		}
//        SimpleAlignment
		
		
	}

	public void simulateSequence(double errorRate, SiteModel siteModel, SubstitutionModel substitutionModel,
			TreeModel treeModel) {

        double substitutionRate = errorRate/(getHaplotypeCount()*getHaplotypeLength()) ;
        System.err.println(substitutionRate);
        double damageRate = 0;
//        SiteModel siteModel = treeLikelihood.getSiteModel();
//        SubstitutionModel substitutionModel = siteModel.getSubstitutionModel();
        
        int[] initialSequence = aMap.getConsensusSequenceState();
        StringBuffer buffer = new StringBuffer();
        for (int i = 0; i < initialSequence.length; i++) {
            buffer.append(Nucleotides.INSTANCE.getChar(    initialSequence[i] ));
        }

        SeqGenExt seqGen = new SeqGenExt(initialSequence, 
                substitutionRate, substitutionModel, siteModel, damageRate);
        
        jebl.evolution.alignments.Alignment jeblAlignment = seqGen.simulate(treeModel);
        List<Sequence> sequenceList = jeblAlignment.getSequenceList();
        for (int j = 0; j < sequenceList.size(); j++) {
			OldHaplotype haplotype = getHaplotype(j);
			haplotype.setSequenceString(sequenceList.get(j).getString());
		}
		
	}
	
	public void swapHaplotypeColumn(int[] posChar){
	
		int[] allOldChars = new int[getHaplotypeCount()];
			
		for (int i = 0; i < getHaplotypeCount(); i++) {
			OldHaplotype haplotype = haplotypes.get(i);
			allOldChars[i] = haplotype.replaceCharAt(posChar[0], posChar[1]);
		}
			
		storeOperationRecord(OldHapOperation.SWAPCOLUMN, posChar, allOldChars);
	
	
	}



	public int[] swapHaplotypeSingleBase(OldHapOperation op, int[] posChar){
	
		swapBaseRecord[POS_INDEX] = posChar[0];
		swapBaseRecord[NEW_CHAR_INDEX] = posChar[1];
	
		swapBaseRecord[HAP_INDEX] = MathUtils.nextInt(getHaplotypeCount());
	
		OldHaplotype haplotype = haplotypes.get(swapBaseRecord[HAP_INDEX]);
		swapBaseRecord[OLD_CHAR_INDEX] = haplotype.replaceCharAt(
				swapBaseRecord[POS_INDEX], swapBaseRecord[NEW_CHAR_INDEX]);
	
		storeOperationRecord(op, swapBaseRecord);
		return swapBaseRecord;
	
	}






	public void swapHaplotypeMultiBases(OldHapOperation op, int[][] allPosChars){
			int hapIndex = MathUtils.nextInt(getHaplotypeCount());
			
			OldHaplotype haplotype = haplotypes.get(hapIndex);
	
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



	



	//	private int swapHaplotypeCharAt(int hapIndex, int[] posChar){
	//		int oldChar = getHaplotype(hapIndex).getChar(posChar[0]); 
	//		getHaplotype(hapIndex).setCharAt(posChar[0], (char) posChar[1]);
	//		return oldChar;
	//		
	//	}
	//	
	//	
	
		public double swapNextDiffBaseFrequency(OldHapOperation op, Parameter frequency) {
			
				checkFrequencyParameter(frequency);
				
		//		int[] swapBaseRecord = new int[4];
				swapBaseRecord[HAP_INDEX] = MathUtils.nextInt(getHaplotypeCount());
				swapBaseRecord[POS_INDEX] = MathUtils.nextInt(getHaplotypeLength());
		
				OldHaplotype haplotype = haplotypes.get(swapBaseRecord[HAP_INDEX]);
				swapBaseRecord[OLD_CHAR_INDEX] = haplotype.getChar(swapBaseRecord[POS_INDEX]);
				
				do{
					swapBaseRecord[NEW_CHAR_INDEX] = setNewCharFromFrequency();
				}while(swapBaseRecord[OLD_CHAR_INDEX]==swapBaseRecord[NEW_CHAR_INDEX]);
		
				haplotype.setCharAt(swapBaseRecord[POS_INDEX], (char) swapBaseRecord[NEW_CHAR_INDEX]);
				
				storeOperationRecord(op, swapBaseRecord);
		//			return swapBaseRecord;
		//			swapInfo.storeOperation(op, swapInfoArray);
		//			int newChar = NUCLEOTIDE_STATES[swapRecord[NEW_CHAR_INDEX]];
		//			int oldChar = NUCLEOTIDE_STATES[swapRecord[OLD_CHAR_INDEX]];
				double logq = getLogqFrequency(swapBaseRecord[OLD_CHAR_INDEX], swapBaseRecord[NEW_CHAR_INDEX]);
				
				return logq;
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
				OldHaplotype haplotype = haplotypes.get(swapRecord[0]);
				haplotype.replaceCharAt(swapRecord[1], swapRecord[3]);
		
		//		swapHaplotypeCharAt(swapArray[0], swapArray[1], swapArray[3]);
			}



		//	public Operation getOperation() {
	//		return swapInfo.getOperation();
	//	}
		public void storeOperationRecord(OldHapOperation op, int[]... opRecord){
			swapInfo.storeOperation(op, opRecord);
		}



	private void storeOperationRecord(OldHapOperation op, int hapIndex,
			int[][] opRecord) {
		swapInfo.storeOperation(op, opRecord[0], opRecord[1]);
		swapInfo.storeHapIndex(hapIndex);
	}
}