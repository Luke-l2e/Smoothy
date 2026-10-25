package org.metagene.genestrip.smoothy.model;

import java.util.Collections;
import java.util.List;

/**
 * What fitting a {@link TopicModel} to one sample gives: the topics with their estimates, and how the
 * fit went.
 */
public class TopicModelResult {
    private final List<Topic> topics;
    private final int iterations;
    private final boolean converged;
    private final long fittedCount;
    private final List<String> topicsWithoutRow;
    private final List<String> uncoveredNodes;

    /**
     * Creates the result.
     *
     * @param topics           the topics, by descending k-mer share
     * @param iterations       the number of EM iterations taken
     * @param converged        whether the fit converged within the allowed iterations
     * @param fittedCount      the total count the fit explains, C of Equation (12)
     * @param topicsWithoutRow admitted nodes that could not be topics because no genome gives them a row
     * @param uncoveredNodes   admitted nodes no topic's lineage reaches, whose counts were left out
     */
    public TopicModelResult(List<Topic> topics, int iterations, boolean converged, long fittedCount,
                            List<String> topicsWithoutRow, List<String> uncoveredNodes) {
        this.topics = Collections.unmodifiableList(topics);
        this.iterations = iterations;
        this.converged = converged;
        this.fittedCount = fittedCount;
        this.topicsWithoutRow = Collections.unmodifiableList(topicsWithoutRow);
        this.uncoveredNodes = Collections.unmodifiableList(uncoveredNodes);
    }

    /**
     * Returns the topics, by descending k-mer share.
     *
     * @return the topics
     */
    public List<Topic> getTopics() {
        return topics;
    }

    /**
     * Returns the number of EM iterations taken.
     *
     * @return the iterations
     */
    public int getIterations() {
        return iterations;
    }

    /**
     * Returns whether the fit converged within the allowed iterations. A fit that did not is still
     * usable - EM raises the likelihood at every step - but slow convergence is a symptom worth
     * reporting: it is what near-identical topics cause.
     *
     * @return whether the fit converged
     */
    public boolean isConverged() {
        return converged;
    }

    /**
     * Returns the total count the fit explains, summed over the admitted nodes some topic reaches.
     *
     * @return C, in k-mers
     */
    public long getFittedCount() {
        return fittedCount;
    }

    /**
     * Returns the admitted nodes that are not topics because no genome contributes a row to them. These
     * are taxa whose genomes all have topics of their own; the paper would give them a row from the
     * genomes the database did not store, which the counts do not include.
     *
     * @return the tax ids of those nodes
     */
    public List<String> getTopicsWithoutRow() {
        return topicsWithoutRow;
    }

    /**
     * Returns the admitted nodes no topic's lineage reaches. Their counts are left out of the fit. With a
     * database and counts that belong together, this list is empty.
     *
     * @return the tax ids of those nodes
     */
    public List<String> getUncoveredNodes() {
        return uncoveredNodes;
    }
}
