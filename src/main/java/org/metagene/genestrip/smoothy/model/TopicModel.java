package org.metagene.genestrip.smoothy.model;

import java.util.*;

/**
 * The topic model of the smoothy paper, Sections 3.1 to 3.6, fitted per sample against the per-genome
 * path counts of one database.
 * <p>
 * It is independent of Genestrip's goals and trees: a genome is given as the tax ids of its path to the
 * root, the node it sits at first, together with its counts &kappa;<sub>h</sub> along that path as
 * {@code KMerGenomePathCountsGoal} measures them. Everything the model needs of the tree follows from
 * those paths. This keeps the arithmetic testable on hand-built trees of a few nodes.
 * <p>
 * A genome the database stored sits at its own genome node. One it did not store, the second and third
 * case of Section 3.2 (p. 21), sits at the taxon it is filed at, or the nearest ancestor the database's
 * tree has: it has no topic of its own and is spoken for by the most specific topic above it, where it
 * adds to the row and to the divisor E<sub>m</sub> like any other member. The model treats both alike;
 * only where the path starts differs.
 * <p>
 * Fitting a sample takes its count c<sub>n</sub> per node and goes through these steps:
 * <ol>
 * <li>Admission, Equation (2): the topics are the nodes whose count reaches the threshold T. Counts at
 * nodes below it are left out of the fit, not folded upwards.</li>
 * <li>Pools, Equation (3): each genome belongs to the most specific topic on its path, G(m) being the
 * genomes of topic m. A genome with no topic on its path belongs to none.</li>
 * <li>Rows, Equation (4): the row of m sums the counts of G(m) over path(m), the nodes from m to the root,
 * normalised to one.</li>
 * <li>Smoothing, Equation (5): each row hands the share &epsilon; to a uniform background over its own
 * path.</li>
 * <li>The fit, Equation (13): expectation-maximisation over the topic shares &theta;, started uniform.</li>
 * <li>Abundances, Equations (14) to (18): each share divided by what one member of the topic is
 * expected to contribute, E<sub>m</sub>, and normalised.</li>
 * </ol>
 * Before the fit, the counts can be attenuated by the k-mer consistency weight of Section 3.3,
 * Equations (7) to (10), see {@link ConsistencyWeight}: it takes effect for the nodes given a
 * {@link NodeConsistency} in {@link #fit(Map, Map)}. Admission still reads the raw counts. The error
 * correction of Equation (6) is not applied, since under maximum likelihood it does not move the
 * estimate.
 * <p>
 * What the counts do not cover, the model cannot either. Given only the stored genomes - the counts
 * measured with {@code smoothy.unstoredGenomes=false} - a taxon whose genomes all have topics of their own
 * gets no row and is reported in {@link TopicModelResult#getTopicsWithoutRow()} rather than fitted.
 */
public class TopicModel {
    /**
     * path(h) per genome h: the tax ids from the node h sits at up to the root.
     */
    private final Map<String, String[]> genomePaths;
    /**
     * κ_h(n) per genome h, aligned with path(h): the k-mer occurrences of h the database holds at each node.
     */
    private final Map<String, long[]> pathCounts;
    /**
     * T, Eq. (2): the count c_n a node needs to be admitted as a topic.
     */
    private final long threshold;
    /**
     * ε, Eq. (5): the share of each row spread uniformly over path(m).
     */
    private final double backgroundShare;
    /**
     * Upper bound on the EM steps of Eq. (13); not in the paper.
     */
    private final int maxIterations;
    /**
     * The fit stops once Σ_m |θ_m(new) - θ_m(old)| falls below this; not in the paper.
     */
    private final double tolerance;
    /**
     * β, τ, γ of Eqs. (8) to (10).
     */
    private final ConsistencyWeight consistencyWeight;

    /**
     * Creates the model for one database, without consistency weighting.
     *
     * @param genomePaths     per genome, the tax ids of its path from its own node to the root
     * @param pathCounts      per genome, its counts &kappa;<sub>h</sub> along that path, aligned with it
     * @param threshold       the admission threshold T of Equation (2), at least one
     * @param backgroundShare the background share &epsilon; of Equation (5), in [0, 1)
     * @param maxIterations   the most EM iterations to take
     * @param tolerance       the total absolute change in the shares below which the fit has converged
     * @throws IllegalArgumentException if a parameter is out of range or a genome's path and counts do
     *                                  not match
     */
    public TopicModel(Map<String, String[]> genomePaths, Map<String, long[]> pathCounts, long threshold,
                      double backgroundShare, int maxIterations, double tolerance) {
        this(genomePaths, pathCounts, threshold, backgroundShare, maxIterations, tolerance, ConsistencyWeight.NONE);
    }

    /**
     * Creates the model for one database.
     *
     * @param genomePaths       per genome, the tax ids of its path from its own node to the root
     * @param pathCounts        per genome, its counts &kappa;<sub>h</sub> along that path, aligned with it
     * @param threshold         the admission threshold T of Equation (2), at least one
     * @param backgroundShare   the background share &epsilon; of Equation (5), in [0, 1)
     * @param maxIterations     the most EM iterations to take
     * @param tolerance         the total absolute change in the shares below which the fit has converged
     * @param consistencyWeight the consistency weighting of Section 3.3, {@link ConsistencyWeight#NONE} for
     *                          none
     * @throws IllegalArgumentException if a parameter is out of range or a genome's path and counts do
     *                                  not match
     */
    public TopicModel(Map<String, String[]> genomePaths, Map<String, long[]> pathCounts, long threshold,
                      double backgroundShare, int maxIterations, double tolerance,
                      ConsistencyWeight consistencyWeight) {
        if (consistencyWeight == null) {
            throw new IllegalArgumentException("The consistency weight must not be null.");
        }
        if (threshold < 1) {
            throw new IllegalArgumentException("The threshold must be at least 1: " + threshold);
        }
        if (backgroundShare < 0 || backgroundShare >= 1) {
            throw new IllegalArgumentException("The background share must be in [0, 1): " + backgroundShare);
        }
        if (maxIterations < 1) {
            throw new IllegalArgumentException("At least one iteration is needed: " + maxIterations);
        }
        for (Map.Entry<String, long[]> genomeCounts : pathCounts.entrySet()) { // genome h -> κ_h along path(h)
            String[] path = genomePaths.get(genomeCounts.getKey());
            if (path == null || path.length != genomeCounts.getValue().length) {
                throw new IllegalArgumentException("Path and counts of genome " + genomeCounts.getKey() + " do not match.");
            }
        }
        this.genomePaths = genomePaths;
        this.pathCounts = pathCounts;
        this.threshold = threshold;
        this.backgroundShare = backgroundShare;
        this.maxIterations = maxIterations;
        this.tolerance = tolerance;
        this.consistencyWeight = consistencyWeight;
    }

    /**
     * Fits the model to one sample, without consistency weighting.
     *
     * @param nodeCounts the sample's count c<sub>n</sub> per node tax id, in k-mers
     * @return the topics and how the fit went
     */
    public TopicModelResult fit(Map<String, Long> nodeCounts) {
        return fit(nodeCounts, Collections.emptyMap());
    }

    /**
     * Fits the model to one sample.
     *
     * @param nodeCounts  the sample's count c<sub>n</sub> per node tax id, in k-mers
     * @param consistency per node tax id, what the consistency weight of Equation (9) needs; a node
     *                    missing here keeps weight one
     * @return the topics and how the fit went
     */
    public TopicModelResult fit(Map<String, Long> nodeCounts, Map<String, NodeConsistency> consistency) {
        // Equation (2): the admitted nodes, in a fixed order so that the fit is reproducible.
        Map<String, Integer> admittedIndex = new LinkedHashMap<>(); // N_s: admitted node tax id -> its index n in the arrays below
        List<String> admittedNodes = new ArrayList<>(); // N_s as a sorted list; the position is the index n
        for (Map.Entry<String, Long> nodeCount : nodeCounts.entrySet()) { // node tax id n -> sample count c_n
            if (nodeCount.getKey() != null && nodeCount.getValue() != null && nodeCount.getValue() >= threshold) {
                admittedNodes.add(nodeCount.getKey());
            }
        }
        admittedNodes.sort(null);
        for (String node : admittedNodes) {
            admittedIndex.put(node, admittedIndex.size());
        }

        // Equation (3): each genome to the most specific admitted node on its path.
        Map<String, List<Member>> pools = new HashMap<>(); // G(m): per topic node m, the genomes it speaks for
        Map<String, String[]> topicPaths = new HashMap<>(); // path(m): per topic node m, the nodes from m up to the root
        List<String> genomes = new ArrayList<>(pathCounts.keySet()); // all genomes h counted, stored or not
        genomes.sort(null);
        for (String genome : genomes) {
            String[] path = genomePaths.get(genome); // path(h)
            for (int pathPosition = 0; pathPosition < path.length; pathPosition++) { // walks up path(h); the first admitted node is h's topic m
                String node = path[pathPosition];
                if (admittedIndex.containsKey(node)) {
                    pools.computeIfAbsent(node, topic -> new ArrayList<>())
                            .add(new Member(pathCounts.get(genome), pathPosition));
                    topicPaths.putIfAbsent(node, Arrays.copyOfRange(path, pathPosition, path.length));
                    break;
                }
            }
        }

        // Equations (4) and (5): one smoothed row per topic, kept only at the admitted nodes, since the
        // likelihood of Equation (1) runs over those alone. Topics are ordered like the admitted nodes.
        List<String> topicNodes = new ArrayList<>(); // the topics m that got a row; the position is the index m
        List<String> topicsWithoutRow = new ArrayList<>(); // admitted nodes with empty G(m) or a row total of 0: not fitted
        List<int[]> rowNodes = new ArrayList<>(); // per topic m, the indices n of N_s ∩ path(m) where its row is kept
        List<double[]> rowValues = new ArrayList<>(); // per topic m, p̃(n | m) = φ_nm of Eq. (5), aligned with rowNodes
        List<double[]> memberHoldsPerTopic = new ArrayList<>(); // per topic m, ρ_h(m) of each member h ∈ G(m), Eq. (14)
        for (String node : admittedNodes) { // node = m
            List<Member> pool = pools.get(node); // G(m)
            String[] path = topicPaths.get(node); // path(m)
            if (pool == null) {
                topicsWithoutRow.add(node);
                continue;
            }
            long[] rowCounts = new long[path.length]; // Σ_{h ∈ G(m)} κ_h(n) per n ∈ path(m): numerator of Eq. (4)
            double[] memberHolds = new double[pool.size()]; // ρ_h(m) = Σ_{n ∈ path(m)} κ_h(n) per member h, Eq. (14)
            for (int memberIndex = 0; memberIndex < pool.size(); memberIndex++) {
                Member member = pool.get(memberIndex);
                for (int pathPosition = 0; pathPosition < rowCounts.length; pathPosition++) {
                    long genomeCount = member.counts[member.offset + pathPosition]; // κ_h(n) at that node n
                    rowCounts[pathPosition] += genomeCount;
                    memberHolds[memberIndex] += genomeCount;
                }
            }
            long rowTotal = 0; // Σ_{n ∈ path(m)} Σ_{h ∈ G(m)} κ_h(n): denominator of Eq. (4)
            for (long rowCount : rowCounts) {
                rowTotal += rowCount;
            }
            if (rowTotal == 0) {
                topicsWithoutRow.add(node);
                continue;
            }
            int admittedOnPath = 0; // |N_s ∩ path(m)|: how many entries of the row the likelihood reads
            for (String pathNode : path) {
                if (admittedIndex.containsKey(pathNode)) {
                    admittedOnPath++;
                }
            }
            int[] keptNodes = new int[admittedOnPath]; // indices n of N_s ∩ path(m)
            double[] keptValues = new double[admittedOnPath];
            int rowEntry = 0;
            for (int pathPosition = 0; pathPosition < path.length; pathPosition++) {
                Integer nodeIndex = admittedIndex.get(path[pathPosition]); // index n of this node of path(m), or null if not in N_s
                if (nodeIndex != null) {
                    keptNodes[rowEntry] = nodeIndex;
                    keptValues[rowEntry] = (1 - backgroundShare) * rowCounts[pathPosition] / rowTotal
                            + backgroundShare / path.length; // Eq. (5): (1 - ε) p(n | m) + ε / |path(m)|
                    rowEntry++;
                }
            }
            topicNodes.add(node);
            rowNodes.add(keptNodes);
            rowValues.add(keptValues);
            memberHoldsPerTopic.add(memberHolds);
        }

        // The counts the fit explains: those at admitted nodes where some topic's row is positive. Lying on a
        // topic's path is not enough: without smoothing the row can be zero there, and so can the mixture
        // (φθ)_n that the EM update divides by.
        boolean[] covered = new boolean[admittedNodes.size()]; // per n: some φ_nm > 0, so (φθ)_n cannot be 0
        for (int m = 0; m < rowNodes.size(); m++) {
            int[] nodes = rowNodes.get(m);
            double[] values = rowValues.get(m);
            for (int rowEntry = 0; rowEntry < nodes.length; rowEntry++) {
                if (values[rowEntry] > 0) {
                    covered[nodes[rowEntry]] = true;
                }
            }
        }
        // Equations (7) to (10): the consistency weight w_n per admitted node, one where it is declined.
        Map<String, NodeWeight> nodeWeights = new LinkedHashMap<>(); // nodes inside the rank gate, by tax id
        double[] weights = new double[admittedNodes.size()]; // w_n per index n
        Arrays.fill(weights, 1);
        if (consistencyWeight.isActive()) {
            for (int n = 0; n < weights.length; n++) {
                String node = admittedNodes.get(n);
                NodeConsistency nodeConsistency = consistency.get(node);
                if (nodeConsistency == null || !nodeConsistency.isGated() || nodeConsistency.getDbKMers() < 1) {
                    continue;
                }
                long count = nodeCounts.get(node);
                double reach = reach(pools.get(node), nodeConsistency.getDbKMers());
                double ratio = ConsistencyWeight.ratio(nodeConsistency.getUniqueKMers(), nodeConsistency.getDbKMers(),
                        count);
                boolean applied = reach >= consistencyWeight.getReachThreshold();
                if (applied) {
                    weights[n] = consistencyWeight.weight(ratio);
                }
                nodeWeights.put(node, new NodeWeight(node, ratio, reach, applied, weights[n]));
            }
        }

        double[] counts = new double[admittedNodes.size()]; // ĉ_n = w_n c_n per index n, Eq. (9); 0 where n is not covered
        long fittedCount = 0; // Σ_n c_n over the covered nodes, unweighted
        double weightedCount = 0; // C = Σ_n ĉ_n over the covered nodes
        List<String> uncoveredNodes = new ArrayList<>(); // admitted nodes no row reaches; their counts leave the fit
        for (int n = 0; n < counts.length; n++) {
            if (covered[n]) {
                long count = nodeCounts.get(admittedNodes.get(n));
                fittedCount += count;
                counts[n] = weights[n] * count;
                weightedCount += counts[n];
            } else {
                uncoveredNodes.add(admittedNodes.get(n));
            }
        }

        // Equation (13), started uniform since the update keeps a zero share at zero.
        int topicCount = topicNodes.size();
        double[] theta = new double[topicCount]; // θ_m = p(m | s): the topic's share of the sample's k-mers
        Arrays.fill(theta, 1.0 / topicCount);
        double[] mixture = new double[counts.length]; // (φθ)_n = Σ_m φ_nm θ_m per node n
        double[] nextTheta = new double[topicCount]; // θ after one EM step
        int iterations = 0;
        // With nothing to fit, the shares stay uniform rather than turning into 0 / 0.
        boolean converged = topicCount <= 1 || weightedCount == 0;
        while (!converged && iterations < maxIterations) {
            iterations++;
            Arrays.fill(mixture, 0);
            for (int m = 0; m < topicCount; m++) {
                int[] nodes = rowNodes.get(m);
                double[] values = rowValues.get(m);
                for (int rowEntry = 0; rowEntry < nodes.length; rowEntry++) {
                    mixture[nodes[rowEntry]] += values[rowEntry] * theta[m];
                }
            }
            double change = 0; // Σ_m |θ_m(new) - θ_m(old)|
            for (int m = 0; m < topicCount; m++) {
                int[] nodes = rowNodes.get(m);
                double[] values = rowValues.get(m);
                double weightedSum = 0; // Σ_n ĉ_n φ_nm / (φθ)_n
                for (int rowEntry = 0; rowEntry < nodes.length; rowEntry++) {
                    int n = nodes[rowEntry];
                    if (counts[n] > 0) {
                        weightedSum += counts[n] * values[rowEntry] / mixture[n];
                    }
                }
                nextTheta[m] = theta[m] * weightedSum / weightedCount; // Eq. (13): θ_m <- θ_m / C · Σ_n ĉ_n φ_nm / (φθ)_n
                change += Math.abs(nextTheta[m] - theta[m]);
            }
            double[] swap = theta;
            theta = nextTheta;
            nextTheta = swap;
            converged = change < tolerance;
        }

        // Equations (14) to (18): the share over what one member is expected to contribute.
        double[] expectedHold = new double[topicCount]; // E_m of Eq. (16): expected k-mers of one member of G(m)
        double abundanceSum = 0; // Σ_m ã_m: denominator of Eq. (18)
        double[] abundance = new double[topicCount]; // ã_m = θ_m / E_m, Eq. (17); a_m once divided by abundanceSum
        for (int m = 0; m < topicCount; m++) {
            double holdSum = 0; // Σ_{h ∈ G(m)} ρ_h(m): the row's total
            double holdSquares = 0; // Σ_{h ∈ G(m)} ρ_h(m)²
            for (double memberHold : memberHoldsPerTopic.get(m)) {
                holdSum += memberHold;
                holdSquares += memberHold * memberHold;
            }
            expectedHold[m] = holdSquares / holdSum; // Eq. (16): Σ_h q_m(h) ρ_h(m) with q_m(h) = ρ_h(m) / holdSum
            abundance[m] = theta[m] / expectedHold[m];
            abundanceSum += abundance[m];
        }

        List<Topic> topics = new ArrayList<>();
        for (int m = 0; m < topicCount; m++) {
            String node = topicNodes.get(m);
            topics.add(new Topic(node, nodeCounts.get(node), pools.get(node).size(), theta[m], expectedHold[m],
                    abundanceSum > 0 ? abundance[m] / abundanceSum : 0, weights[admittedIndex.get(node)]));
        }
        topics.sort(Comparator.comparingDouble(Topic::getKmerShare).reversed());
        return new TopicModelResult(topics, iterations, converged, fittedCount, weightedCount, topicsWithoutRow,
                uncoveredNodes, new ArrayList<>(nodeWeights.values()));
    }

    /**
     * The reach f<sub>n</sub> of Equation (8), read off the occurrence counts as Appendix F describes: the
     * least share &kappa;<sub>h</sub>(n) / u<sub>n,d</sub> of the node's stored k-mers a member of its pool
     * holds at the node itself, over the members whose pool weight q<sub>n</sub>(h) of Equation (16) is at
     * least &gamma;.
     *
     * @param pool    G(n), the genomes the node speaks for, or null if none
     * @param dbKMers u<sub>n,d</sub>, at least one
     * @return f<sub>n</sub>, or 0 if the pool is empty or holds nothing, which declines the weight
     */
    private double reach(List<Member> pool, long dbKMers) {
        if (pool == null) {
            return 0;
        }
        double[] holds = new double[pool.size()]; // ρ_h(n) = Σ of κ_h over path(n), Eq. (14)
        double holdSum = 0;
        for (int i = 0; i < holds.length; i++) {
            Member member = pool.get(i);
            for (int k = member.offset; k < member.counts.length; k++) {
                holds[i] += member.counts[k];
            }
            holdSum += holds[i];
        }
        if (holdSum == 0) {
            return 0;
        }
        double reach = Double.POSITIVE_INFINITY;
        for (int i = 0; i < holds.length; i++) {
            if (holds[i] / holdSum >= consistencyWeight.getMemberFloor()) { // q_n(h) ≥ γ
                Member member = pool.get(i);
                reach = Math.min(reach, ((double) member.counts[member.offset]) / dbKMers); // κ_h(n) / u_n,d
            }
        }
        // With γ above every q_n(h) no member is left; the test is then declined rather than passed.
        return reach == Double.POSITIVE_INFINITY ? 0 : reach;
    }

    /**
     * A genome as a member of a topic's pool: its counts and where the topic's node sits on its path.
     */
    private static class Member {
        /**
         * κ_h along path(h), index 0 being h's own node.
         */
        private final long[] counts;
        /**
         * The position of the topic m on path(h), so counts[offset + i] = κ_h(i-th node of path(m)).
         */
        private final int offset;

        private Member(long[] counts, int offset) {
            this.counts = counts;
            this.offset = offset;
        }
    }
}
