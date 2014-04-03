package srp.haplotypes;



public class HaplotypeOperationRecord {

	/*
	 * Record how/which move/operation is performed
	*/
	public static final int SWAPBASE_HAP_INDEX = 0;
	public static final int SWAPBASE_POS_INDEX = 1;
	public static final int SWAPBASE_NEW_CHAR_INDEX = 2;
	public static final int SWAPBASE_OLD_CHAR_INDEX = 3;
	
//	private Operations data = new Operations();
//	private int hapIndex;
	
	private HaplotypeOperation operation;
	private int spectrumIndex;
	private int singleIndex;
	private double[] delta;
	private int[] allSiteIndexs;
	private int[] recombinationPositionIndex;
	private int[] recombinationSpectrumIndex;
	private int[] allSpectrumIndexs;

	

		
	
	public HaplotypeOperationRecord() {
		operation = HaplotypeOperation.NONE;

	}


	public HaplotypeOperation getOperation(){
//		operation = SpectrumOperation.NONE;
		return operation;
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


	public int getSpectrumIndex() {
		return spectrumIndex;
	}


	public int getSingleIndex() {
		return singleIndex;
	}

	public int[] getAllSiteIndexs() {
		return allSiteIndexs;
	}


	public double[] getDelta() {
		// for unit test only
		return delta;
	}


	public void setOperation(HaplotypeOperation operation) {
		this.operation = operation;
		
	}

	@Deprecated
	public void setRecord(HaplotypeOperation op, int spectrumIndex, int siteIndex, double[] delta) {
		//Single
		setOperation(op);
		this.spectrumIndex = spectrumIndex;
		this.singleIndex = siteIndex;
		this.delta = delta;
		
	}


	public void setRecord(HaplotypeOperation op, int columnIndex, double[] delta) {
		//Column
		setOperation(op);
		this.singleIndex = columnIndex;
		this.delta = delta;
				
	}
//	public void setRecord(SpectrumOperation op, int spectrumIndex, int[] siteIndexs) {
//		//Swap
//		setOperation(op);
//		this.spectrumIndex = spectrumIndex;
//		this.allSiteIndexs = siteIndexs;
//	}


	public void setRecord(HaplotypeOperation op, int spectrumIndex, int[] siteIndexs, double[] delta) {
		//Multi //Swap
		setOperation(op);
		this.spectrumIndex = spectrumIndex;
		this.allSiteIndexs = siteIndexs;
		this.delta = delta;
	}
	public void setRecord(HaplotypeOperation op, int spectrumIndex, int siteIndex) {
		//Single
		setOperation(op);
		this.spectrumIndex = spectrumIndex;
		this.singleIndex= siteIndex;
//		this.delta = delta;
	}


	public void setRecord(HaplotypeOperation op, int[] recombinationSpectrumIndex,
			int[] recombinationPositionIndex) {
		setOperation(op);
		this.recombinationSpectrumIndex = recombinationSpectrumIndex;
		this.recombinationPositionIndex = recombinationPositionIndex;

	}


	public void setRecord(HaplotypeOperation op, int[] hapIndexs, int siteIndex) {
		setOperation(op);
		this.singleIndex = siteIndex;
		this.allSpectrumIndexs = hapIndexs;

	}



	public int[] getRecombinationPositionIndex() {
		return recombinationPositionIndex;
	}


	public int[] getRecombinationSpectrumIndex() {
		return recombinationSpectrumIndex;
	}


	public int[] getAllSpectrumIndexs() {
		return allSpectrumIndexs;
	}


	
}
