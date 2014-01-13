package srp.spectrum;


public enum SpectrumOperation{
	NONE(0),PASS(0),
	
	SINGLE_DELTA(1), COLUMN_DELTA(2), MULTI_DELTA(3), RECOMBINATION(4), 
	
////	UNIFORMSWAPBASE(6),
//	SWAPMULTI(3),
//	 
//	SWAPCOLUMN(2), 
//	
//	SWAPSECTION(5), SWITCH_PROB(6),  
//	
//	TREE(10)
	;
	
	private int code;

	private SpectrumOperation(int i){
		code = i;
	}
	public int getCode(){
		return code;
	}
}
