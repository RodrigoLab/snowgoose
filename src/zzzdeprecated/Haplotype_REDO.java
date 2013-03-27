package zzzdeprecated;



import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;
import dr.util.Attributable;
import dr.util.Identifiable;

import java.util.Iterator;
@Deprecated
public class Haplotype_REDO implements Identifiable, Attributable {

	    /**
	     * Empty constructor.
	     */
	    public Haplotype_REDO() {
	        haplotypeString = new StringBuilder();
	    }

	    /**
	     * Constructor with initial haplotype string.
	     *
	     * @param haplotype a string representing the haplotype
	     */
	    public Haplotype_REDO(String haplotype) {
	        haplotypeString = new StringBuilder();
	        setHaplotypeString(haplotype);
	    }

	    /**
	     * Clone constructor
	     *
	     * @param haplotype_REDO the haplotype to clone
	     */
	    public Haplotype_REDO(Haplotype_REDO haplotype_REDO) {
	        // should clone taxon as well!
	        this(haplotype_REDO.getTaxon(), haplotype_REDO.getHaplotypeString());
	    }

	    /**
	     * Constructor with taxon and haplotype string.
	     *
	     * @param taxon    the haplotype's taxon
	     * @param haplotype the haplotype's symbol string
	     */
	    public Haplotype_REDO(Taxon taxon, String haplotype) {
	        haplotypeString = new StringBuilder();
	        setTaxon(taxon);
	        setHaplotypeString(haplotype);
	    }



	    /**
	     * @return the length of the haplotype.
	     */
	    public int getLength() {
	        return haplotypeString.length();
	    }

	    /**
	     * @return a String containing the haplotype.
	     */
	    public String getHaplotypeString() {
	        return haplotypeString.toString();
	    }

	    /**
	     * @return a char containing the state at index.
	     */
	    public char getCharAt(int index) {
	        return haplotypeString.charAt(index);
	    }
	    
	    public final void setCharAt(int index, char ch) {

	        haplotypeString.setCharAt(index, ch);
	    }
	    
	    /**
	     * @return the state at site index.
	     */
	    public int getState(int index) {
	        return dataType.getState(haplotypeString.charAt(index));
	    }

	    /**
	     */
	    public final void setState(int index, int state) {

	        haplotypeString.setCharAt(index, dataType.getChar(state));
	    }

	    /**
	     * Characters are copied from the haplotype into the destination character array dst.
	     */
	    public void getChars(int srcBegin, int srcEnd, char[] dst, int dstBegin) {
	        haplotypeString.getChars(srcBegin, srcEnd, dst, dstBegin);
	    }

	    /**
	     * Set the DataType of the haplotype.
	     */
	    public void setDataType(DataType dataType) {
	        this.dataType = dataType;
	    }
	    
	    /**
	     * @return the DataType of the haplotype.
	     */
	    public DataType getDataType() {
	        return dataType;
	    }
	    
	    /**
	     * Set the DataType of the haplotype.
	     */
	    public DataType guessDataType() {
	        return DataType.guessDataType(haplotypeString.toString());
	    }

	    /**
	     * Set the haplotype using a string.
	     */
	    public void setHaplotypeString(String haplotype) {
	        haplotypeString.setLength(0);
	        haplotypeString.append(haplotype.toUpperCase());
	    }

	    /**
	     * Append a string to the haplotype.
	     */
	    public void appendHaplotypeString(String haplotype) {
	        haplotypeString.append(haplotype);
	    }

	    /**
	     * Insert a string into the haplotype.
	     */
	    public void insertHaplotypeString(int offset, String haplotype) {
	        haplotypeString.insert(offset, haplotype);
	    }

	    /**
	     * Sets a taxon for this haplotype.
	     *
	     * @param taxon the taxon.
	     */
	    public void setTaxon(Taxon taxon) {
	        this.taxon = taxon;
	    }

	    /**
	     * @return the taxon for this haplotype.
	     */
	    public Taxon getTaxon() {
	        return taxon;
	    }

	    // **************************************************************
	    // Attributable IMPLEMENTATION
	    // **************************************************************

	    private Attributable.AttributeHelper attributes = null;

	    /**
	     * Sets an named attribute for this object.
	     *
	     * @param name  the name of the attribute.
	     * @param value the new value of the attribute.
	     */
	    public void setAttribute(String name, Object value) {
	        if (attributes == null)
	            attributes = new Attributable.AttributeHelper();
	        attributes.setAttribute(name, value);
	    }

	    /**
	     * @param name the name of the attribute of interest.
	     * @return an object representing the named attributed for this object.
	     */
	    public Object getAttribute(String name) {
	        if (attributes == null)
	            return null;
	        else
	            return attributes.getAttribute(name);
	    }

	    /**
	     * @return an iterator of the attributes that this object has.
	     */
	    public Iterator<String> getAttributeNames() {
	        if (attributes == null)
	            return null;
	        else
	            return attributes.getAttributeNames();
	    }

	    // **************************************************************
	    // Identifiable IMPLEMENTATION
	    // **************************************************************

	    protected String id = null;

	    /**
	     * @return the id.
	     */
	    public String getId() {
	        return id;
	    }

	    /**
	     * Sets the id.
	     */
	    public void setId(String id) {
	        this.id = id;
	    }

	    // **************************************************************
	    // INSTANCE VARIABLES
	    // **************************************************************

	    protected Taxon taxon = null;
	    protected StringBuilder haplotypeString = null;
	    protected DataType dataType = null;
	}





