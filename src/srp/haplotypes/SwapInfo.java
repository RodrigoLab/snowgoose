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
	
	int[] swapBase = new int[4];
	
	private Operation operation;
	private Deque<int[]> swapMulti;

	
	private int[] swapHaplotypeRecord;
	
	public SwapInfo() {
		operation = Operation.NONE;
//		Operation a = Operation.SWAPBASE;
//		Operation swapbase = Operation.SWAPBASE;
//		Oper
//		Operation.valueOf(arg0)
	}
	
	public Operation getOperation(){
		return operation;
	}

	public void storeOperation(Operation op, Object swapRecord){
		this.operation = op;
		int[] tempIntArray;
		switch (operation) {
			case NONE:
				break;
			case SWAPBASE:
				swapBase = (int[]) swapRecord;
				break;
			case UNIFORMSWAPBASE:
				swapBase = (int[]) swapRecord;
				break;
//			case SWAPCOLUMN:
//				
//				break;
			case SWAPMULTI:
				if (swapRecord==null){
					swapMulti = new ArrayDeque<int[]>(); 
//					swapMulti = new ArrayList<int[]>();
				}
				else{
//					tempIntArray = new int[4];
					tempIntArray = (int[]) swapRecord;
//					for (int i = 0; i < tempIntArray.length; i++) {
//						tempIntArray[i] = (Integer) swapInfo[i];
//					}
					swapMulti.add( tempIntArray);
				}
				break;
			case SWAPSECTION:
				swapHaplotypeRecord = (int[]) swapRecord;
				break;
			case RECOMB:
				swapHaplotypeRecord = (int[]) swapRecord;
				break;
			default:
				throw new IllegalArgumentException("Unknown operation type: "+op);
			}
		
		}

	public int[] getSwapInfoSWAPBASE(){
		return swapBase;
	
	}

	public Deque<int[]> getSwapInfoSWAPMULTI(){
		return swapMulti;
	}


	public class InvalidOperationException extends Exception {

		public InvalidOperationException(String message) {
			super(message);
		}
	}

	public int[] getSwapHaplotypeRecord() {
		
		return swapHaplotypeRecord;
	}	
}
	

