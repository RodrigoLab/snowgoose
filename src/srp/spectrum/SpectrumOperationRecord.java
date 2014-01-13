package srp.spectrum;

import java.util.Arrays;


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
	private int columnIndex;
	private double[] delta;
	private int[] allSiteIndexs;
	private int[] recombinationPositionIndex;
	private int[] recombinationSpectrumIndex;

	
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


	public int getSpectrumIndex() {
		return spectrumIndex;
	}
//
//
//	public void setSpectrumIndex(int spectrumIndex) {
//		this.spectrumIndex = spectrumIndex;
//	}


	public int getColumnIndex() {
		return columnIndex;
	}

	public int[] getAllSiteIndexs() {
	
		return allSiteIndexs;
	}


	public double[] getDelta() {
		return delta;
	}


	public void setOperation(SpectrumOperation operation) {
		this.operation = operation;
		
	}

	@Deprecated
	public void setRecord(SpectrumOperation op, int spectrumIndex, int siteIndex, double[] delta) {
		//Single
		setOperation(op);
		this.spectrumIndex = spectrumIndex;
		this.columnIndex = siteIndex;
		this.delta = delta;
		
	}


	public void setRecord(SpectrumOperation op, int columnIndex, double[] delta) {
		//Column
		setOperation(op);
		this.columnIndex = columnIndex;
		this.delta = delta;
				
	}


	public void setRecord(SpectrumOperation op, int spectrumIndex, int[] siteIndexs, double[] delta) {
		//Multi
		setOperation(op);
		this.spectrumIndex = spectrumIndex;
		this.allSiteIndexs = siteIndexs;
		this.delta = delta;
		
		
	}



	public void setRecord(SpectrumOperation op, int[] recombinationSpectrumIndex,
			int[] recombinationPositionIndex) {
		setOperation(op);
		this.recombinationSpectrumIndex = recombinationSpectrumIndex;
		this.recombinationPositionIndex = recombinationPositionIndex;

	}


	public int[] getRecombinationPositionIndex() {
		return recombinationPositionIndex;
	}


	public int[] getRecombinationSpectrumIndex() {
		return recombinationSpectrumIndex;
	}

	
}
