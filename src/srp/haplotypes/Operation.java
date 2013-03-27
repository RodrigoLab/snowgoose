package srp.haplotypes;

public enum Operation {
	NONE(0),
	SWAPBASE(1), SWAPCOLUMN(2);
	
	private int code;
	private Operation(int i){
		code = i;
	}
	public int getCode(){
		return code;
	}
	

}
