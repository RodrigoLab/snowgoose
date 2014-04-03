package srp.haplotypes;


public enum HaplotypeOperation{
	NONE(0),PASS(0), 
	FULL(1),
	
	SINGLE(1),
	MULTI(2),
	
	
//	RECOMBINATION(4), 
//	COLUMN(5),
//	DIRICHLET(6), 
//	
	;
	
	private int type;

	private HaplotypeOperation(int i){
		type = i;
	}
	public int getCode(){
		return type;
	}
}
