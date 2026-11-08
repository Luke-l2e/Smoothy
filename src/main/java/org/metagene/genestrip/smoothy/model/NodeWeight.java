package org.metagene.genestrip.smoothy.model;

/**
 * How the consistency weight of Section 3.3 came out for one admitted node inside the rank gate: the
 * quantities of Equations (7) to (9) and the weight the node's count entered the fit with.
 * <p>
 * Kept for every gated node, also where the weight was declined, since the distribution of the reach
 * and of the ratio over a database is what &tau;, &gamma; and &beta; are to be fixed against (the
 * {@code todo}s at the end of Section 3.3).
 */
public class NodeWeight {
    private final String node;
    private final double ratio;
    private final double reach;
    private final boolean applied;
    private final double weight;

    /**
     * Creates the node's weight record.
     *
     * @param node    the node's tax id
     * @param ratio   the consistency ratio r<sub>n,s</sub> of Equation (7)
     * @param reach   the reach f<sub>n</sub> of Equation (8)
     * @param applied whether f<sub>n</sub> &ge; &tau;, i.e. whether the weight was applied
     * @param weight  the weight w<sub>n</sub> the count entered the fit with, one if not applied
     */
    public NodeWeight(String node, double ratio, double reach, boolean applied, double weight) {
        this.node = node;
        this.ratio = ratio;
        this.reach = reach;
        this.applied = applied;
        this.weight = weight;
    }

    /**
     * Returns the node's tax id.
     *
     * @return the tax id
     */
    public String getNode() {
        return node;
    }

    /**
     * Returns the consistency ratio of Equation (7).
     *
     * @return r<sub>n,s</sub>
     */
    public double getRatio() {
        return ratio;
    }

    /**
     * Returns the reach of Equation (8).
     *
     * @return f<sub>n</sub>
     */
    public double getReach() {
        return reach;
    }

    /**
     * Returns whether the weight was applied, the reach having passed &tau;.
     *
     * @return whether f<sub>n</sub> &ge; &tau;
     */
    public boolean isApplied() {
        return applied;
    }

    /**
     * Returns the weight the node's count entered the fit with.
     *
     * @return w<sub>n</sub>
     */
    public double getWeight() {
        return weight;
    }

    @Override
    public String toString() {
        return "NodeWeight[" + node + ", ratio=" + ratio + ", reach=" + reach + ", applied=" + applied + ", weight="
                + weight + "]";
    }
}
