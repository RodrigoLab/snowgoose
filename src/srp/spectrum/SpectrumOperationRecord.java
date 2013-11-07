package srp.spectrum;

import srp.haplotypes.Operation;

public class SpectrumOperationRecord {

	/*
	 * Record how/which move/operation is performed
	*/
	public static final int SWAPBASE_HAP_INDEX = 0;
	public static final int SWAPBASE_POS_INDEX = 1;
	public static final int SWAPBASE_NEW_CHAR_INDEX = 2;
	public static final int SWAPBASE_OLD_CHAR_INDEX = 3;
	
//	private Operations data = new Operations();
//	private int hapIndex;
	
	private int[] swapBaseRecord;// hapIndex, posIndex, newChar, oldChar
	
	private int[] swapHaplotypeSectoinRecord;
	private int[][] allPosChars;
	private int[] allOldChars;
	private int[] posChars;
	private SpectrumOperation operation;
	private int spectrumIndex;
	private int siteIndex;
	private double[] delta;

	
	/*
		public void storeOperation(Operation op, int[] swapRecord){
			
			this.operation = op;
	
			switch (operation) {
				case NONE:
					break;
				case SWAPSINGLE:
					swapBase = swapRecord;
					break;
	
	//			case SWAPMULTI:
	//				hapIndex = (Integer) swapRecord[0];
	//				allPosChars = (int[][]) swapRecord[1];
	//
	//				break;
				case SWAPSECTION:
					swapHaplotypeRecord = swapRecord;
					break;
	//			case SWAPCOLUMN:
	//				posChars = (int[]) swapRecord[0];
	//				allOldChars = (int[]) swapRecord[1];
				default:
					throw new IllegalArgumentException("Unknown operation type: "+op);
				}
			
		}
		*/
		
	
	public SpectrumOperationRecord() {
		// TODO Auto-generated constructor stub
		operation = SpectrumOperation.NONE;
	}


	public SpectrumOperation getOperation(){
//		operation = SpectrumOperation.NONE;
		return operation;
	}
	/*
	public void storeOperation(Operation op, int[] swapRecord){
		
		this.operation = op;

		switch (operation) {
			case NONE:
				break;
			case SWAPSINGLE:
				swapBase = swapRecord;
				break;

//			case SWAPMULTI:
//				hapIndex = (Integer) swapRecord[0];
//				allPosChars = (int[][]) swapRecord[1];
//
//				break;
			case SWAPSECTION:
				swapHaplotypeRecord = swapRecord;
				break;
//			case SWAPCOLUMN:
//				posChars = (int[]) swapRecord[0];
//				allOldChars = (int[]) swapRecord[1];
			default:
				throw new IllegalArgumentException("Unknown operation type: "+op);
			}
		
	}
	*/
	
	
	public void storeOperation(SpectrumOperation op, int[]... swapRecord){
		
		operation = op;

		switch (operation) {
			case NONE:
				break;
//			case SWAPSINGLE:
//				swapBaseRecord = swapRecord[0];
//				break;
//
//			case SWAPMULTI:
//				allPosChars[0] = swapRecord[0];
//				allPosChars[1] = swapRecord[1];
//
//			case SWAPSECTION:
//				swapHaplotypeSectoinRecord = swapRecord[0];
//				break;
//			case SWAPCOLUMN:
//				posChars = swapRecord[0];
//				allOldChars = swapRecord[1];
//				break;
			case PASS:
				break;
			default:
				throw new IllegalArgumentException("Unknown operation type: "+op);
			}
		
	}
	
//	public void storeHapIndex(int hapIndex2) {
//		hapIndex = hapIndex2;
//	
//		
//	}

	public int[] getSwapInfoSWAPBASE(){
		return swapBaseRecord;
	
	}

	public int[][] getSwapInfoSWAPMULTI(){
		return allPosChars;
	}

	
	public int[][] getSwapInfoSWAPCOLUMN() {
		int[][] record = new int[][]{posChars,	allOldChars}; 
		return record;
	}

	public int[] getSwapInfoSWAPSECTION() {
		
		return swapHaplotypeSectoinRecord;
	}

//	public int getHapIndex(){
//		return hapIndex;
//	}

	
	public class InvalidOperationException extends Exception {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public InvalidOperationException(String message) {
			super(message);
		}
	}


	public void setRecord(SpectrumOperation operation, int spectrumIndex, int siteIndex, double[] delta) {
		this.operation = operation;
		this.spectrumIndex = spectrumIndex;
		this.siteIndex = siteIndex;
		this.delta = delta;
		
	}


	public int getSpectrumIndex() {
		return spectrumIndex;
	}
//
//
//	public void setSpectrumIndex(int spectrumIndex) {
//		this.spectrumIndex = spectrumIndex;
//	}


	public int getSiteIndex() {
		return siteIndex;
	}
//
//
//	public void setSiteIndex(int siteIndex) {
//		this.siteIndex = siteIndex;
//	}


	public double[] getDelta() {
		return delta;
	}


	public void setOperation(SpectrumOperation operation) {
		this.operation = operation;
		
	}


	
}
