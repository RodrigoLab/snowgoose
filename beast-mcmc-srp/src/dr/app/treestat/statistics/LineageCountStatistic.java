package dr.app.treestat.statistics;

import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;

/**
 * Returns the total time in the genealogy in which exactly k lineages are present.
 *
 * @version $Id: IntervalKStatistic.java,v 1.2 2005/09/28 13:50:56 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class LineageCountStatistic extends AbstractTreeSummaryStatistic {

	public LineageCountStatistic() {
		this.t = 1.0;
	}

    public void setDouble(double value) {
        this.t = value;
    }

	public double[] getSummaryStatistic(Tree tree) {

		TreeIntervals intervals = new TreeIntervals(tree);

		double totalTime = 0.0;
		for (int i = 0; i < intervals.getIntervalCount(); i++) {
			totalTime += intervals.getInterval(i);
			if (totalTime > t) {
				return new double[] { intervals.getLineageCount(i) };
			}
		}
		return new double[] { 1.0 };
	}

	public String getSummaryStatisticName() {
        return "LineageCount(" + t + ")";
    }

	public String getSummaryStatisticDescription() {
        return getSummaryStatisticName() + " is the number of lineages that exists in the genealogy at " +
            "time " + t + ".";
    }

	public String getSummaryStatisticReference() { return FACTORY.getSummaryStatisticReference(); }
	public boolean allowsPolytomies() { return FACTORY.allowsPolytomies(); }
	public boolean allowsNonultrametricTrees() { return FACTORY.allowsNonultrametricTrees(); }
	public boolean allowsUnrootedTrees() { return FACTORY.allowsUnrootedTrees(); }
	public Category getCategory() { return FACTORY.getCategory(); }

	public static final Factory FACTORY = new Factory() {

		public TreeSummaryStatistic createStatistic() {
			return new LineageCountStatistic();
		}

		public String getSummaryStatisticName() {
			return "LineageCount(t)";
		}

		public String getSummaryStatisticDescription() {
			return getSummaryStatisticName() + " is the number of lineages that exists in the genealogy at " +
				"time t.";
		}

		public String getSummaryStatisticReference() {
			return "-";
		}

		public String getValueName() { return "The time (t):"; }
		public boolean allowsPolytomies() { return true; }

		public boolean allowsNonultrametricTrees() { return true; }

		public boolean allowsUnrootedTrees() { return false; }

		public Category getCategory() { return Category.POPULATION_GENETIC; }

        public boolean allowsWholeTree() { return true; }

        public boolean allowsCharacter() { return false; }

        public boolean allowsCharacterState() { return false; }

        public boolean allowsTaxonList() { return false; }

        public boolean allowsInteger() { return false; }

        public boolean allowsDouble() { return true; }
    };

	double t = 1.0;
}
