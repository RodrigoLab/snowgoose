package srp.operator.spectrum;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import com.google.common.primitives.Ints;

import srp.evolution.OperationType;
import srp.evolution.shortreads.AlignmentMapping;
import srp.evolution.spectrum.SpectraParameter;
import srp.evolution.spectrum.SpectrumAlignmentModel;
import srp.evolution.spectrum.SpectrumLogger;
import dr.inference.model.Bounds;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

public abstract class AbstractSpectrumOperator extends AbstractCoercableOperator {

	public static final int DIMENSION = SpectraParameter.DIMENSION;
    public static final Bounds<Double> BOUNDS = SpectraParameter.SPECTRA_BOUNDS;
    public static final double BOUNDS_LOWER = BOUNDS.getLowerLimit(0);
    public static final double BOUNDS_UPPER = BOUNDS.getUpperLimit(0);
    

    public OperationType OP;
    
    protected static final int MIN_BASE = 1;
    protected SpectrumAlignmentModel spectrumModel;

    protected final int[] fixSpectrumIndexArray;
    
	protected final int spectrumCount;
	protected final int spectrumLength;
	
	
	
	private int[] spectrumLengthArray;


//	protected int swapLength;
//
//	protected int[][] allPosChars;
//
//	private double autoOptimize;
//	private double scaleFactor;
	
	public AbstractSpectrumOperator(SpectrumAlignmentModel spectrumModel, CoercionMode mode) {
		super(mode);
		this.spectrumModel = spectrumModel;
		spectrumCount = this.spectrumModel.getSpectrumCount();
		spectrumLength = this.spectrumModel.getSpectrumLength();
		spectrumLengthArray = new int[spectrumLength];
		for (int i = 0; i < spectrumLengthArray.length; i++) {
			spectrumLengthArray[i] = i;
		}
		fixSpectrumIndexArray = new int[spectrumCount];
        for (int i = 0; i < fixSpectrumIndexArray.length; i++) {
        	fixSpectrumIndexArray[i] = i;
		}
	}


	public abstract OperationType getSpectrumOperation();

	HashSet<Integer> generated = new HashSet<Integer>();
	public long time = 0;
	public long time2 = 0;
	public long time3 = 0;
	private long timeStart = 0;

	public int[] generateUniqueSites(int m) {
		int[] siteIndexs;
//		timeStart = System.nanoTime();
//		randomSiteHashSet(m);
//		time += (System.nanoTime() - timeStart);
		
//		timeStart  = System.nanoTime();		
		siteIndexs = randomSampleSites(m, spectrumLengthArray);
//		time2 += (System.nanoTime() - timeStart);
		
//		timeStart  = System.nanoTime();		
//		siteIndexs = randomSiteFloyd(m);
//		time3 += (System.nanoTime() - timeStart);
		
		return siteIndexs;
	}
	public int[] randomSiteHashSet(int m, int total){
		generated.clear();
		while (generated.size() < m) //time: 1.2/1.3
		{
		    Integer next = MathUtils.nextInt(total);
		    generated.add(next);
		}
		
//		int[] siteIndexs = new int[noSample]; //time: ~0.4/0.5
//		int count = 0;
//		for (Integer i : generated) {
//			siteIndexs[count] = i;
//			count++;
//		}
		int[] siteIndexs = Ints.toArray(generated); //time: 0.6
		return siteIndexs;
	}
	
//	public int[] randomSampleSites(int m){ //time:0.5
//		int[] sites = new int[m];
//	    for(int i=0;i<m;i++){
//	        int pos = i + MathUtils.nextInt(spectrumLength - i);
//	        int tmp = spectrumLengthArray[pos];
//	        spectrumLengthArray[pos] = spectrumLengthArray[i];
//	        spectrumLengthArray[i] = tmp;
//	        sites[i] = tmp;
//	        
//	    }
//	    return sites;
//	}
	
	public static int[] randomSampleSites(int m, int[] sampleArray){ //time:0.5
		int[] sites = new int[m];
		int length = sampleArray.length;
		for (int i = 0; i < m; i++) {
	        int pos = i + MathUtils.nextInt(length - i);
//	        T tmp = items.get(pos);
//	        items.set(pos, items.get(i));
//	        items.set(i, tmp);
	        int tmp = sampleArray[pos];
	        sampleArray[pos] = sampleArray[i];
	        sampleArray[i] = tmp;
	        sites[i] = tmp;
	        
	    }
//	    return items.subList(0, m);
	    return sites;
	}

	public int[] randomSiteFloyd(int m){
	    generated.clear();
	    int n = spectrumLength;
	    for(int i=n-m;i<n;i++){
	        int pos = MathUtils.nextInt(i+1);
	        int item = spectrumLengthArray[pos];
	        if (generated.contains(item))
	        	generated.add(spectrumLengthArray[i]);
	        else
	        	generated.add(item);
	    }
	    int[] siteIndexs = Ints.toArray(generated); 
	    return siteIndexs;
	}
	
	public static <T> List<T> randomSampleSwap(List<T> items, int m){
	    for(int i=0;i<m;i++){
	        int pos = i + MathUtils.nextInt(items.size() - i);
	        T tmp = items.get(pos);
	        items.set(pos, items.get(i));
	        items.set(i, tmp);
	    }
	    return items.subList(0, m);
	}
	
	public static <T> Set<T> randomSampleFloyd(List<T> items, int m){
	    HashSet<T> res = new HashSet<T>(m);
	    int n = items.size();
	    for(int i=n-m;i<n;i++){
	        int pos = MathUtils.nextInt(i+1);
	        T item = items.get(pos);
	        if (res.contains(item))
	            res.add(items.get(i));
	        else
	            res.add(item);
	    }
	    return res;
	}
	
	
	
	protected int getAnotherDimension(int dim1) {
		// get any two dims and swap
        int dim2;// = dim1;
        do {
            dim2 = MathUtils.nextInt(DIMENSION);
        }while (dim1 == dim2);
		return dim2;
	}


	
//	public AbstractMultiBasesOperator(HaplotypeModel haplotypeModel, int swapLength, CoercionMode mode) {
//		super(mode);
//		
//		this.haplotypeModel = haplotypeModel;
//		this.swapLength = swapLength;
//		haplotypeLength = this.haplotypeModel.getHaplotypeLength();
//		alignmentMapping = this.haplotypeModel.getAlignmentMapping();
//		
//		allPosChars = new int[2][haplotypeLength];
//	
//		
//		checkParameterIsValid();
//		
//		scaleFactor = (int) (haplotypeLength*0.01);
//
//		if (scaleFactor <1) {
//			scaleFactor = 1;
//		}
//		
//		convertToAutoOptimize(this.swapLength);
//		
//	}

//	@Override
//	public double getCoercableParameter() {
//	    return autoOptimize;
//	}
//
//	@Override
//	public void setCoercableParameter(double autoOpt) {
//		convertFromAutoOptimizeToValue(autoOpt);
//	    
//	}

//	private void convertFromAutoOptimizeToValue(double autoOpt) {
//	    	autoOptimize = autoOpt;
//			swapLength =  1 + (int) Math.exp(autoOptimize*scaleFactor);
////			System.out.println(autoOptimize +"\t"+ Math.exp(autoOptimize*scaleFactor));
//			
////			System.out.print("A=" + swapLength + "\t" + autoOptimize + "\t" +
////					"accept: " + getAcceptCount()/(double)getCount() + "\t"  );
//			
//			checkParameterIsValid();
//			
////			System.out.print("newL:"+swapLength+" ");
//	//		System.out.print("A\t" + swapFragmentLength + "\t" + autoOptimize + "\t"  );
//	    }
//
//	private double convertToAutoOptimize(int length) {
//		swapLength = length;
//		checkParameterIsValid();
//		autoOptimize = Math.log(swapLength - 1)/scaleFactor;
//	    return autoOptimize;
//	}
//
//	private void checkParameterIsValid() {
//		if (swapLength > spectrumLength){
//			swapLength = spectrumLength;
//		}
//	}
//	
//	protected void resetAllPosChars() {
//		for (int i = 0; i < allPosChars.length; i++) {
//			Arrays.fill(allPosChars[i], -1);
//		}
//	}
//
//
//	@Override
//	public double getRawParameter() {
//		return swapLength;
//	}

//
//	@Override
//	public String getPerformanceSuggestion() {
//		return "";
//	}


}



