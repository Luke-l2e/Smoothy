package org.metagene.genestrip.smoothy.model;

/**
 * The k-mer consistency weight of Section 3.3 of the smoothy paper, Equations (7) to (10): a node whose
 * counts show fewer distinct k-mers than uniform sampling from the database's k-mers there would give
 * enters the fit with fewer counts.
 * <p>
 * The weight of node n is
 * <pre>
 *   w_n = w(r_n,s)  if n is at rank &rho; or below and f_n &ge; &tau;,
 *   w_n = 1         otherwise,                                        (Equation 9)
 * </pre>
 * with
 * <ul>
 * <li>r<sub>n,s</sub> = u<sub>n,s</sub> / E<sub>n,s</sub><br>
 * E<sub>n,s</sub> = u<sub>n,d</sub> (1 - (1 - 1/u<sub>n,d</sub>)<sup>c<sub>n</sub></sup>), the
 * consistency ratio of Equation (7), u<sub>n,s</sub> being the distinct k-mers the sample shows at n
 * and u<sub>n,d</sub> the distinct k-mers the database stores there;</li>
 * <li>f<sub>n</sub> the reach of Equation (8): the least share of n's stored k-mers any member of n's
 * pool holds, over the members whose weight q<sub>n</sub>(h) is at least &gamma;;</li>
 * <li>w(r) = max(0, 1 + &beta; ln min(r, 1)), Equation (10).</li>
 * </ul>
 * The rank gate &rho; is applied by the caller, which knows the tree; see {@link NodeConsistency}.
 * <p>
 * Unlike the paper, the reach is read off the occurrence counts &kappa;<sub>h</sub>(n) rather than off
 * HyperLogLog sketches of the distinct k-mers - the variant of Appendix F. It is exact where a genome
 * holds each k-mer once and is inflated by the copy number of repeated k-mers otherwise, so it may admit
 * a node to the weighting that a sketch would decline.
 */
public class ConsistencyWeight {
    /**
     * The weighting that never attenuates: every node keeps weight one.
     */
    public static final ConsistencyWeight NONE = new ConsistencyWeight(0, 0.9, 0.01);

    private final double severity; // β, Eq. (10): how hard an inconsistent node is attenuated; 0 = off
    private final double reachThreshold; // τ, Eq. (9): the reach f_n a node needs to be weighted at all
    private final double memberFloor; // γ, Eq. (8): the pool weight q_n(h) a genome needs to count as member

    /**
     * Creates the weighting.
     *
     * @param severity       &beta; of Equation (10), at least zero; zero turns the weighting off
     * @param reachThreshold &tau; of Equation (9), in [0, 1]
     * @param memberFloor    &gamma; of Equation (8), in [0, 1]
     * @throws IllegalArgumentException if a parameter is out of range
     */
    public ConsistencyWeight(double severity, double reachThreshold, double memberFloor) {
        if (!(severity >= 0) || Double.isInfinite(severity)) {
            throw new IllegalArgumentException("The severity must be finite and at least 0: " + severity);
        }
        if (!(reachThreshold >= 0 && reachThreshold <= 1)) {
            throw new IllegalArgumentException("The reach threshold must be in [0, 1]: " + reachThreshold);
        }
        if (!(memberFloor >= 0 && memberFloor <= 1)) {
            throw new IllegalArgumentException("The member floor must be in [0, 1]: " + memberFloor);
        }
        this.severity = severity;
        this.reachThreshold = reachThreshold;
        this.memberFloor = memberFloor;
    }

    /**
     * Returns whether the weighting can attenuate a node at all.
     *
     * @return whether &beta; is above zero
     */
    public boolean isActive() {
        return severity > 0;
    }

    /**
     * Returns &beta; of Equation (10).
     *
     * @return the severity
     */
    public double getSeverity() {
        return severity;
    }

    /**
     * Returns &tau; of Equation (9).
     *
     * @return the reach threshold
     */
    public double getReachThreshold() {
        return reachThreshold;
    }

    /**
     * Returns &gamma; of Equation (8).
     *
     * @return the member floor
     */
    public double getMemberFloor() {
        return memberFloor;
    }

    /**
     * The distinct k-mers expected at a node under uniform sampling, E<sub>n,s</sub> of Equation (7).
     *
     * @param dbKMers the distinct k-mers u<sub>n,d</sub> the database stores at the node, at least one
     * @param count   the k-mers c<sub>n</sub> the sample reports there
     * @return E<sub>n,s</sub>
     */
    public static double expectedUniqueKMers(long dbKMers, long count) {
        // u (1 - (1 - 1/u)^c), written with expm1 and log1p so that it stays exact for large u.
        return -dbKMers * Math.expm1(count * Math.log1p(-1d / dbKMers));
    }

    /**
     * The consistency ratio r<sub>n,s</sub> of Equation (7).
     *
     * @param uniqueKMers the distinct k-mers u<sub>n,s</sub> the sample shows at the node
     * @param dbKMers     the distinct k-mers u<sub>n,d</sub> the database stores there, at least one
     * @param count       the k-mers c<sub>n</sub> the sample reports there, at least one
     * @return r<sub>n,s</sub>
     */
    public static double ratio(long uniqueKMers, long dbKMers, long count) {
        return uniqueKMers / expectedUniqueKMers(dbKMers, count);
    }

    /**
     * The weight w(r) of Equation (10).
     *
     * @param ratio the consistency ratio r
     * @return max(0, 1 + &beta; ln min(r, 1))
     */
    public double weight(double ratio) {
        if (!(ratio < 1)) {
            return 1; // Also catches NaN: no evidence against the node.
        }
        return Math.max(0, 1 + severity * Math.log(ratio));
    }
}
