package srp.dr.evolution.alignment;

import srp.dr.evolution.datatype.ShortReads;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;

@Deprecated
public class ShortReadFragmentsAlignment extends SimpleAlignment {

	private final static DataType DATA_TYPE= ShortReads.INSTANCE;
	
    
	public ShortReadFragmentsAlignment(){
		super();
		setDataType(DATA_TYPE);
	}
}
