package srp.evolution;



public class OperationRecord {

	/*
	 * Record how/which move/operation is performed
	*/
	
//	private Operations data = new Operations();
//	private int hapIndex;
	
	private OperationType operation;
	private int spectrumIndex;
	private int singleIndex;
	
	private int[] allSiteIndexs;
	private int[] recombinationPositionIndex;
	private int[] recombinationSpectrumIndex;
	private int[] allSpectrumIndexs;

	@Deprecated private double[] delta;

		
	
	public OperationRecord() {
		operation = OperationType.NONE;

	}


	public OperationType getOperation(){
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


	public void setOperation(OperationType operation) {
		this.operation = operation;
		
	}

	@Deprecated
	public void setRecord(OperationType op, int spectrumIndex, int siteIndex, double[] delta) {
		//Single
		setOperation(op);
		this.spectrumIndex = spectrumIndex;
		this.singleIndex = siteIndex;
		this.delta = delta;
		
	}

	@Deprecated
	public void setRecord(OperationType op, int columnIndex, double[] delta) {
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

	@Deprecated
	public void setRecord(OperationType op, int spectrumIndex, int[] siteIndexs, double[] delta) {
		//Multi //Swap
		setOperation(op);
		this.spectrumIndex = spectrumIndex;
		this.allSiteIndexs = siteIndexs;
		this.delta = delta;
	}
	public void setRecord(OperationType op, int spectrumIndex, int[] siteIndexs) {
		//Multi //Swap
		setOperation(op);
		this.spectrumIndex = spectrumIndex;
		this.allSiteIndexs = siteIndexs;
	}
	public void setRecord(OperationType op, int spectrumIndex, int siteIndex) {
		//Single
		setOperation(op);
		this.spectrumIndex = spectrumIndex;
		this.singleIndex= siteIndex;
//		this.delta = delta;
	}


	public void setRecord(OperationType op, int[] recombinationSpectrumIndex,
			int[] recombinationPositionIndex) {
		setOperation(op);
		this.recombinationSpectrumIndex = recombinationSpectrumIndex;
		this.recombinationPositionIndex = recombinationPositionIndex;

	}


	public void setRecord(OperationType op, int[] hapIndexs, int siteIndex) {
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
