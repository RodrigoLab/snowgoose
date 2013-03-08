package srp.haplotypes;

public enum Operation {
	SWAPBASE(0), SWAPCOLUMN(1);
	
	private int code;
	private Operation(int i){
		code = i;
	}
	public int getCode(){
		return code;
	}
	

}
