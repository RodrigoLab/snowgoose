package srp.haplotypes;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.LinkedList;
public class SwapInfo {

	/*
	 * Record how/which move/operation is performed
	*/

	private Operation operation;
	private int hapIndex;
	
	private int[] swapBase;// hapIndex, posIndex, newChar, oldChar
	private int[] swapHaplotypeSectoinRecord;
	private int[][] allPosChars;
	private int[] allOldChars;
	private int[] posChars;
	
	public SwapInfo() {
		operation = Operation.NONE;

	}
	
	public Operation getOperation(){
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
	
	
	public void storeOperation(Operation op, int[]... swapRecord){
		
		this.operation = op;

		switch (operation) {
			case NONE:
				break;
			case SWAPSINGLE:
				swapBase = swapRecord[0];
				break;

//			case SWAPMULTI:
//				hapIndex = (Integer) swapRecord[0];
//				allPosChars = (int[][]) swapRecord[1];
//
//				break;
			case SWAPSECTION:
				swapHaplotypeSectoinRecord = swapRecord[0];
				break;
			case SWAPCOLUMN:
				posChars = swapRecord[0];
				allOldChars = swapRecord[1];
				break;
			default:
				throw new IllegalArgumentException("Unknown operation type: "+op);
			}
		
	}
	
	public void storeOperation(Operation op, Object... swapRecord){
		this.operation = op;

		switch (operation) {
			case NONE:
				break;
//			case SWAPSINGLE:
//				swapBase = (int[]) swapRecord[0];
//				break;

			case SWAPMULTI:
				hapIndex = (Integer) swapRecord[0];
				allPosChars = (int[][]) swapRecord[1];

				break;
//			case SWAPSECTION:
//				swapHaplotypeRecord = (int[]) swapRecord[0];
//				break;
//			case SWAPCOLUMN:
//				posChars = (int[]) swapRecord[0];
//				allOldChars = (int[]) swapRecord[1];
			default:
				throw new IllegalArgumentException("Unknown operation type: "+op);
			}
		
	}

	public int[] getSwapInfoSWAPBASE(){
		return swapBase;
	
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

	public int getHapIndex(){
		return hapIndex;
	}


	public class InvalidOperationException extends Exception {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public InvalidOperationException(String message) {
			super(message);
		}
	}	
}
	

