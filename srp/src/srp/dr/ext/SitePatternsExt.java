package srp.dr.ext;

import java.util.Arrays;

import srp.evolution.OperationRecord;
import srp.haplotypes.HaplotypeModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SiteList;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.util.TaxonList;


public class SitePatternsExt extends SitePatterns {

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
	
	public SitePatternsExt(SiteList siteList) {
		this(siteList, 0, siteList.getSiteCount() -1, 1, true, true);
	}
	
	public SitePatternsExt(SiteList siteList, int from, int to, int every,
			boolean strip, boolean unique) {
		super(siteList, from, to, every, strip, unique);

      if (this.from <= -1)
          this.from = 0;

      if (this.to <= -1)
          this.to = siteList.getSiteCount() - 1;

      if (this.every <= 0)
          this.every = 1;

	}

//	public SitePatternsExt(HaplotypeModel haplotypeModel, TaxonList taxa, int from,
//			int to, int every, boolean strip) {
////		Alignment alignment = haplotypes.getAlignment();
//		this(haplotypeModel.getAlignment(), taxa, from, to, every, strip);
//	}
	public void updateAlignment(HaplotypeModel haplotypeModel){

      this.siteList = haplotypeModel;

//      patternCount = 0;

//      invariantCount = 0;
      int[] pattern;

//      int site = 0;
		OperationRecord record = haplotypeModel.getOperationRecord();
		int hapIndex = record.getSpectrumIndex();
		int[] siteIndex = record.getAllSiteIndexs();
//System.out.println(hapIndex);
		switch (record.getOperation()) {
		case SINGLE:

			break;
		case MULTI:
//			 for (int site = from; site <= to; site += every) {
			for (int site : siteIndex) {
				int[] oldPattern = haplotypeModel.getStoredSitePattern(site);
//System.out.println(patternCount +"\t"+ Arrays.toString(pattern) +"\t"+ sitePatternIndices[site]);				
				
//System.out.println(patternCount);
				pattern = haplotypeModel.getSitePattern(site);
//				System.out.println(site +"\t"+  Arrays.toString(oldPattern) +"\t"+  Arrays.toString(pattern));
				removePatternExt(oldPattern);
				sitePatternIndices[site] = addPatternExt(pattern);
//System.out.println(patternCount +"\t"+ Arrays.toString(pattern) +"\t"+ sitePatternIndices[site]);
//	System.out.println();
//	System.exit(-1);
			}
//			System.exit(-1);
			break;
		case COLUMN:

			break;
		case RECOMBINATION:

			break;

		default:
			throw new IllegalArgumentException("Invalid operation type "
					+ record.getOperation());

		}
	}

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
					System.out.println("weigths at "+ i +"==" +"\t"+ weights[i] +"\t What to do now? patternCount--??" );
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
System.out.println(" new pattern. "+patternCount);
		int index = patternCount;
		patterns[index] = pattern;
		weights[index] = 1.0;
		patternCount++;

		return index;
	}

	//
//	
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
	 * 
	 */
	private static final long serialVersionUID = 8376785905234747979L;

//	SitePatterns patterns2 = new SitePatterns(alignment, null, 0, -1, 1, true);
	
}
