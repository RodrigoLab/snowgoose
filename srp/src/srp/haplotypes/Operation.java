package srp.haplotypes;


public enum Operation{
	NONE(0),PASS(0),
	
	SWAPSINGLE(1), 
	
//	UNIFORMSWAPBASE(6),
	SWAPMULTI(3),
	 
	SWAPCOLUMN(2), 
	
	SWAPSECTION(5), SWITCH_PROB(6),  
	
	TREE(10);
	
	private int code;

	private Operation(int i){
		code = i;
	}
	public int getCode(){
		return code;
	}
}
