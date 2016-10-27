package srp.dr.ext;

import java.util.Arrays;

import srp.evolution.OperationRecord;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SiteList;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.util.TaxonList;



/**
 * 
 * @author sw167
 *
 *	There are several possible implementation
 *	with remove/add patterns
 *	 - need to keep track of the number index and pattern count
 *	 - need to prune patterns
 * 1346966067	134.69660670000002	sitePatternExt.updateAlignment(haplotypeModel);
 * 980860917	98.08609170000001	updatePatternListExt();
 * 
 *	with specturmLike idea, number of patterns == number of site. No pooling things together
 *   - might be slow for long (very long) sequences with repeat patterns
 *   - very fast updating alignment
 *   - basically 
 *   	Arrays.fill(weights, 1.0);
        patternCount = siteCount;
 *   
 *   78226249	7.8226249			sitePatternExt.updateAlignment(haplotypeModel);
 *   1163143736	116.31437360000001	updatePatternListExt();
 *   
 *   update only the change
 *   ~7-10 sitePatternExt.updateAlignment(haplotypeModel);
 *   ~7-10 updatePatternListExt();
 */
public class SitePatternsExt extends SitePatterns {

	private static final long serialVersionUID = 8376785905234747979L;
	private static final boolean UNIQUE_FALSE = false;
	private static final boolean STRIP_FALSE = false;
	
	@Deprecated
	public SitePatternsExt(Alignment alignment, TaxonList taxa, int from,
			int to, int every, boolean strip) {
		super(alignment, taxa, from, to, every, strip);

      if (this.from <= -1)
          this.from = 0;

      if (this.to <= -1)
          this.to = siteList.getSiteCount() - 1;

      if (this.every <= 0)
          this.every = 1;

	}

	public SitePatternsExt(SiteList siteList, int from, int to, int every){
		super(siteList, from, to, every, STRIP_FALSE, UNIQUE_FALSE);
		System.err.println("skiping with \"every\" is not tested, don't think it's going to work");
	}

	
	public SitePatternsExt(SiteList siteList) {
		super(siteList, 0, siteList.getSiteCount() - 1, 1, STRIP_FALSE, UNIQUE_FALSE);
	}
	

	public void updateAlignment(HaplotypeModel haplotypeModel) {

		this.siteList = haplotypeModel;

		OperationRecord record = haplotypeModel.getOperationRecord();
		int hapIndex = record.getSpectrumIndex();
//		int[] pattern;
		int state;
		Haplotype haplotype;
		
		switch (record.getOperation()) {
		case SINGLE:
			haplotype = haplotypeModel.getHaplotype(hapIndex);
			int site = record.getSingleIndex();
			state = haplotype.getState(site);
			patterns[site][hapIndex] = state;

			break;
		case MULTI:
			haplotype = haplotypeModel.getHaplotype(hapIndex);
			int[] siteIndex = record.getAllSiteIndexs();
			for (int s : siteIndex) {
				state = haplotype.getState(s);
				patterns[s][hapIndex] = state;
			}
			break;
		case COLUMN:
			site = record.getSingleIndex();
			for ( hapIndex = 0; hapIndex < haplotypeModel.getHaplotypeCount(); hapIndex++) {
				haplotype = haplotypeModel.getHaplotype(hapIndex);
				state = haplotype.getState(site);
				patterns[site][hapIndex] = state;
			}
			break;
		case RECOMBINATION:
			int[] twoPositions = record.getRecombinationPositionIndex();
			int[] twoHapIndex = record.getRecombinationSpectrumIndex();
			haplotype = haplotypeModel.getHaplotype(twoHapIndex[0]);
			Haplotype haplotype1 = haplotypeModel.getHaplotype(twoHapIndex[1]);
			for (int s = twoPositions[0]; s < twoPositions[1]; s++) {
				state = haplotype.getState(s);
				patterns[s][twoHapIndex[0]] = state;
				state = haplotype1.getState(s);
				patterns[s][twoHapIndex[1]] = state;
			}
			
			break;
		case FULL:
//			System.out.println("Update sitepatternExt FULL");
			updateAlignment(haplotypeModel, 0);
			break;
		default:
			throw new IllegalArgumentException("Invalid operation type "
					+ record.getOperation());

			// case MULTI:
			// for (int site : siteIndex) {
			// pattern = haplotypeModel.getStoredSitePattern(site);
			// removePatternExt(pattern);
			// pattern = haplotypeModel.getSitePattern(site);
			// sitePatternIndices[site] = addPatternExt(pattern);
			// }
		}
	}

	@Deprecated
	private void setAllPatterns() {

        if (siteList == null) {
            return;
        }

        if (from <= -1)
            from = 0;

        if (to <= -1)
            to = siteList.getSiteCount() - 1;

        if (every <= 0)
            every = 1;

//        siteCount = ((to - from) / every) + 1;
        patternCount = 0;
        invariantCount = 0;

        sitePatternIndices = new int[siteCount];
        
//        System.out.println(siteList.getPatternLength());
//        System.out.println(siteList.getTaxonCount());
        int scaler = siteList.getPatternCount();
//        scaler=
        int maxPatternCount = siteCount ;//* scaler;
        pruningThreshold = siteCount -10;
        patterns = new int[maxPatternCount][];
        weights = new double[maxPatternCount];
        Arrays.fill(weights, 1.0);
        patternCount = siteCount;
        pruningThreshold = siteCount -10;
        int site = 0;

        for (int i = from; i <= to; i += every) {
            int[] pattern = siteList.getSitePattern(i);
            patterns[i] = pattern;
            sitePatternIndices[site] = i;
//            if (!strip || !isInvariant(pattern) ||
//                    (!isGapped(pattern) &&
//                            !isAmbiguous(pattern) &&
//                            !isUnknown(pattern))) {
//
//                sitePatternIndices[site] = addPattern(pattern);
//
//            }  else {
//                sitePatternIndices[site] = -1;
//            }
            site++;
        }
    }

	@Deprecated
	private void removePatternExt(int[] pattern) {
//		System.out.print("Remove match at: ");
		for (int i = 0; i < patternCount; i++) {
			if (comparePatterns(patterns[i], pattern)) {
//				System.out.println("pattern "+i);
				weights[i] -= 1.0;
				if(weights[i] <= 0){
//					weights[i] = weights[patternCount-1];
//					for (int j = 0; j < patterns[i].length; j++) {
//						patterns[i][j] = patterns[patternCount-1][j];
//					}
//					patternCount--;
//					System.out.println("weigths at "+ i +"==" +"\t"+ weights[i] +"\t What to do now? patternCount--??" );
					prune = true;
				}
				break;
			}
		}
//		System.out.println("cant find match??");
//		System.out.println(Arrays.toString(pattern));
//		for (int i = 0; i < patternCount; i++) {
//			System.out.println(Arrays.toString(patterns[i]));
//		}
//		System.out.println(patterns[patternCount+1]);
//		System.exit(-1);
		if (isInvariant(pattern)) {
			invariantCount--;
		}
	}

	@Deprecated
	private int addPatternExt(int[] pattern) {
//		System.out.print("add pattern: ");
		for (int i = 0; i < patternCount; i++) {
			if (unique && comparePatterns(patterns[i], pattern)) {
//				System.out.println(" at i "+i);
				weights[i] += 1.0;
				return i;
			}
		}

		if (isInvariant(pattern)) {
			invariantCount++;
		}
//System.out.println(" new pattern. "+patternCount);
		int index = patternCount;
		patterns[index] = pattern;
		weights[index] = 1.0;
		patternCount++;
//        System.out.println(patternCount);
//      if(patternCount == patterns.length){
//      	prunePatterns();
//      }


		return index;
	}

	@Deprecated
    private void prunePatterns() {
    	
    	int pruneCount = 0;
    	double[] tempWeights = new double[weights.length];
    	
//    	int[] tempPatterns= new int[patternCount];
    	System.arraycopy(weights, 0, tempWeights, 0, weights.length);
//    	System.arraycopy(patterns, 0, tempPatterns, 0, patternCount);
    	for (int i = 0; i < patternCount; i++) {

    		if(tempWeights[i] > 0.0 ){
    			
    			weights[pruneCount] = tempWeights[i];
    			patterns[pruneCount] = patterns[i];
    			pruneCount++;
    		}
            
        }
    	if(pruneCount<patternCount){
        	System.out.println("pattern conut "+ patternCount +"\tPruneCount"+ pruneCount +"\tdelta: "+ (patternCount-pruneCount) );
    	}
    	patternCount = pruneCount;
    	
//    	if()
//    	System.out.println("pattern conut "+ patternCount);
//    	System.out.println("new pruneCount "+pruneCount);
		
	}

	//
//	
    @Deprecated
	public void updateAlignment(Alignment HaplotypeModel, int x){

//        setPatterns(alignment, from, to, every);

        this.siteList = HaplotypeModel;
//        this.from = from;
//        this.to = to;
//        this.every = every;

//        if (siteList == null) {
//            return;
//        }

//        if (from <= -1)
//            from = 0;

//        if (to <= -1)
//            to = siteList.getSiteCount() - 1;

//        if (every <= 0)
//            every = 1;

//        siteCount = ((to - from) / every) + 1;

        patternCount = 0;

//        patterns = new int[siteCount][];

//        sitePatternIndices = new int[siteCount];
//        weights = new double[siteCount];

        invariantCount = 0;
        int[] pattern;

        int site = 0;


        for (int i = from; i <= to; i += every) {
            pattern = siteList.getSitePattern(i);

            if (!strip || !isInvariant(pattern) ||
                    (!isGapped(pattern) &&
                            !isAmbiguous(pattern) &&
                            !isUnknown(pattern))) {

                sitePatternIndices[site] = addPattern(pattern);

            }  else {
              sitePatternIndices[site] = -1;
            }
            site++;
        }
	}
    @Deprecated
    private int addPattern(int[] pattern) {
	
	    for (int i = 0; i < patternCount; i++) {
	
	        if (unique && comparePatterns(patterns[i], pattern)) {
	
	            weights[i] += 1.0;
	            return i;
	        }
	    }
	
	    if (isInvariant(pattern)) {
	        invariantCount++;
	    }
	
	    int index = patternCount;
	    patterns[index] = pattern;
	    weights[index] = 1.0;
	    patternCount++;
	
	    return index;
	}

	/**
     * @return true if the pattern is invariant
     */
    private boolean isGapped(int[] pattern) {
        int len = pattern.length;

        for (int i = 0; i < len; i++) {
            if (getDataType().isGapState(pattern[i])) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return true if the pattern is invariant
     */
    private boolean isAmbiguous(int[] pattern) {
        int len = pattern.length;

        for (int i = 0; i < len; i++) {
            if (getDataType().isAmbiguousState(pattern[i])) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return true if the pattern is invariant
     */
    private boolean isUnknown(int[] pattern) {
        int len = pattern.length;

        for (int i = 0; i < len; i++) {
            if (getDataType().isUnknownState(pattern[i])) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return true if the pattern is invariant
     */
    private boolean isInvariant(int[] pattern) {
        int len = pattern.length;

        int state = pattern[0];
        for (int i = 1; i < len; i++) {
            if (pattern[i] != state) {
                return false;
            }
        }

        return true;
    }
	private boolean prune;
	private int pruningThreshold;

	public void restoreState(HaplotypeModel haplotypeModel) {
		updateAlignment(haplotypeModel);
//		OperationRecord record = haplotypeModel.getOperationRecord();
//		int hapIndex = record.getSpectrumIndex();
//		int[] siteIndex = record.getAllSiteIndexs();
//		int[] pattern;
//		switch (record.getOperation()) {
//		case SINGLE:
//		case COLUMN:
//			int site = record.getSingleIndex();
//			pattern = haplotypeModel.getStoredSitePattern(site);
//			patterns[site] = pattern;
//
//			break;
//		
//		}
	}
//	SitePatterns patterns2 = new SitePatterns(alignment, null, 0, -1, 1, true);

	public void storeState() {
		// TODO Auto-generated method stub
		
	}
	
}
