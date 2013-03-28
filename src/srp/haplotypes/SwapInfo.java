package srp.haplotypes;


public class SwapInfo {

	/*
	 * Record how/which move/operation is performed
	*/
	
	int[] swapBase = new int[4];
	private Operation operation;
	
	public SwapInfo() {
		operation = Operation.NONE;
//		Operation a = Operation.SWAPBASE;
//		Operation swapbase = Operation.SWAPBASE;
//		Oper
//		Operation.valueOf(arg0)
	}
	
	public void storeOperation(Operation op, int[] swapInfoOld){
		this.operation = op;
		switch (op) {
			case SWAPBASE:
				for (int i = 0; i < swapBase.length; i++) {
					swapBase[i] = (Integer) swapInfoOld[i];
				}
				break;
			case SWAPCOLUMN:
				
				break;
			default:
				break;
			}
		
		}
//	
	public int[] getSwapInfoIntArray(){//FIXME
		switch (operation) {
		case SWAPBASE:
			return swapBase;
			
		case SWAPCOLUMN:
			
			break;
		default:
			break;
		}
		return swapBase;
	
	}
	
	public Operation getOperation(){
		return operation;
	}
	
}
	

