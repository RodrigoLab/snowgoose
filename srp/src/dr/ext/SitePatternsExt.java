package dr.ext;

import dr.evolution.alignment.Alignment;
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
//	public SitePatternsExt(HaplotypeModel haplotypeModel, TaxonList taxa, int from,
//			int to, int every, boolean strip) {
////		Alignment alignment = haplotypes.getAlignment();
//		this(haplotypeModel.getAlignment(), taxa, from, to, every, strip);
//	}
	public void updateAlignment(Alignment HaplotypeModel){

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
        weights = new double[siteCount];

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
