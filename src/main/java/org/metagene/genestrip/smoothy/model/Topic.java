package org.metagene.genestrip.smoothy.model;

/**
 * One topic of a fitted {@link TopicModel}: a node admitted for the sample, with what the fit estimated
 * for it.
 * <p>
 * Two shares are kept apart because they answer different questions (Section 3.6 of the smoothy
 * paper). The k-mer share &theta;<sub>m</sub> is what Equation (1) estimates: the share of the sample's
 * k-mers the topic accounts for, and so a share of sequence. The abundance a<sub>m</sub> divides it by
 * what one organism of the topic contributes, E<sub>m</sub> of Equation (16), and so is a share of
 * organisms; a genome twice as long has twice the k-mer share at equal abundance.
 */
public class Topic {
    private final String node;
    private final long count;
    private final int genomes;
    private final double kmerShare;
    private final double expectedHold;
    private final double abundance;

    /**
     * Creates the topic.
     *
     * @param node         the tax id of the topic's node
     * @param count        the count c<sub>m</sub> the sample reports at the node
     * @param genomes      the number of genomes the topic speaks for, |G(m)| of Equation (3)
     * @param kmerShare    the estimated share &theta;<sub>m</sub> of the sample's k-mers
     * @param expectedHold the expected hold E<sub>m</sub> of one member, in k-mers
     * @param abundance    the relative abundance a<sub>m</sub> of Equation (18)
     */
    public Topic(String node, long count, int genomes, double kmerShare, double expectedHold, double abundance) {
        this.node = node;
        this.count = count;
        this.genomes = genomes;
        this.kmerShare = kmerShare;
        this.expectedHold = expectedHold;
        this.abundance = abundance;
    }

    /**
     * Returns the tax id of the topic's node.
     *
     * @return the node's tax id
     */
    public String getNode() {
        return node;
    }

    /**
     * Returns the count the sample reports at the topic's node.
     *
     * @return c<sub>m</sub>, in k-mers
     */
    public long getCount() {
        return count;
    }

    /**
     * Returns how many genomes the topic speaks for: one at a genome node, the genomes without a topic of
     * their own beneath a taxon.
     *
     * @return |G(m)|
     */
    public int getGenomes() {
        return genomes;
    }

    /**
     * Returns the estimated share of the sample's k-mers this topic accounts for.
     *
     * @return &theta;<sub>m</sub>
     */
    public double getKmerShare() {
        return kmerShare;
    }

    /**
     * Returns how many k-mers one organism of this topic is expected to contribute.
     *
     * @return E<sub>m</sub>, in k-mers
     */
    public double getExpectedHold() {
        return expectedHold;
    }

    /**
     * Returns the relative abundance of the topic: its share of organisms rather than of sequence.
     *
     * @return a<sub>m</sub>
     */
    public double getAbundance() {
        return abundance;
    }

    @Override
    public String toString() {
        return "Topic[" + node + ", count=" + count + ", genomes=" + genomes + ", kmerShare=" + kmerShare
                + ", abundance=" + abundance + "]";
    }
}
