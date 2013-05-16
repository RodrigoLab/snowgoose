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
	
	private int[] swapBase;// hapIndex, posIndex, newChar, oldChar

	private int[] swapHaplotypeRecord;

	private int hapIndex;
	private int[][] allPosChars;
	
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

	public void storeOperation(Operation op, Object... swapRecord){
		this.operation = op;

		switch (operation) {
			case NONE:
				break;
			case SWAPSINGLE:
				swapBase = (int[]) swapRecord[0];
				break;

			case SWAPMULTI:
				hapIndex = (Integer) swapRecord[0];
				allPosChars = (int[][]) swapRecord[1];

//				int[] allNewChars = (int[]) swapRecord[2];
				
//				System.out.println(Arrays.toString(allNewChars));
//				if (swapRecord==null){
//					swapMulti = new ArrayDeque<int[]>(); 
////					swapMulti = new ArrayList<int[]>();
//				}
//				else{
////					tempIntArray = new int[4];
//					tempIntArray = (int[]) swapRecord[0];
////					for (int i = 0; i < tempIntArray.length; i++) {
////						tempIntArray[i] = (Integer) swapInfo[i];
////					}
//					swapMulti.add( tempIntArray);
//				}
				break;
			case SWAPSECTION:
				swapHaplotypeRecord = (int[]) swapRecord[0];
				break;

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

	
	public int getHapIndex(){
		return hapIndex;
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
	

