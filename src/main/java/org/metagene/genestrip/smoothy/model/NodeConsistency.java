package org.metagene.genestrip.smoothy.model;

/**
 * What the consistency weight of Section 3.3 needs to know about one node for one sample, beside its
 * count c<sub>n</sub>.
 */
public class NodeConsistency {
    private final long uniqueKMers;
    private final long dbKMers;
    private final boolean gated;

    /**
     * Creates the node's consistency data.
     *
     * @param uniqueKMers the distinct k-mers u<sub>n,s</sub> the sample shows at the node
     * @param dbKMers     the distinct k-mers u<sub>n,d</sub> the database stores at the node
     * @param gated       whether the node is at rank &rho; or below, i.e. inside the rank gate of
     *                    Equation (9); outside it the node keeps weight one
     */
    public NodeConsistency(long uniqueKMers, long dbKMers, boolean gated) {
        this.uniqueKMers = uniqueKMers;
        this.dbKMers = dbKMers;
        this.gated = gated;
    }

    /**
     * Returns the distinct k-mers the sample shows at the node.
     *
     * @return u<sub>n,s</sub>
     */
    public long getUniqueKMers() {
        return uniqueKMers;
    }

    /**
     * Returns the distinct k-mers the database stores at the node.
     *
     * @return u<sub>n,d</sub>
     */
    public long getDbKMers() {
        return dbKMers;
    }

    /**
     * Returns whether the node is inside the rank gate of Equation (9).
     *
     * @return whether the node is at rank &rho; or below
     */
    public boolean isGated() {
        return gated;
    }
}
