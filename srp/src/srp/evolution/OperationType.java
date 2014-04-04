package srp.evolution;


public enum OperationType{
	NONE(0),PASS(0), 
	FULL(1),
	
	SINGLE(1),  MULTI(2),
	
	COLUMN(3),
	
	RECOMBINATION(4), 
	
	DIRICHLET(6), 
	
	SWAP_SUBCOLUMN(7), 
//	MULTI(7), SWAP_COLUMN(7), SWAP_SUBCOLUMN(8),
	
//	DELTA_SINGLE(2), DELTA_COLUMN(2), DELTA_MULTI(2),
//	
//	RECOMBINATION(4), 
//	
//	DIRICHLET(6), 
//	SWAP_SINGLE(7), SWAP_MULTI(7), SWAP_COLUMN(7), SWAP_SUBCOLUMN(8),
//	
////	UNIFORMSWAPBASE(6),
//	SWAPMULTI(3),
//	 
//	SWAPCOLUMN(2), 
//	
//	SWAPSECTION(5), SWITCH_PROB(6),  
//	
//	TREE(10)
	;
	
	private int type;

	private OperationType(int i){
		type = i;
	}
	public int getCode(){
		return type;
	}
}
