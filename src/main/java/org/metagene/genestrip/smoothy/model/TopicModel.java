package org.metagene.genestrip.smoothy.model;

import java.util.*;

/**
 * The topic model of the smoothy paper, Sections 3.1 to 3.6, fitted per sample against the per-genome
 * path counts of one database.
 * <p>
 * It is independent of Genestrip's goals and trees: a genome is given as the tax ids of its path to the
 * root, its own node first, together with its counts &kappa;<sub>h</sub> along that path as
 * {@code KMerGenomePathCountsGoal} measures them. Everything the model needs of the tree follows from
 * those paths. This keeps the arithmetic testable on hand-built trees of a few nodes.
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
 * The correction and consistency weighting of the counts (Section 3.3) are not applied: under maximum
 * likelihood the error correction of Equation (6) does not move the estimate, and the consistency weight
 * is optional (Table 5).
 * <p>
 * What the counts do not cover, the model cannot either: the paper also pools genomes the database did
 * not store into the topics at the taxa, while the counts hold only stored genomes. A taxon whose
 * genomes all have topics of their own therefore gets no row and is reported in
 * {@link TopicModelResult#getTopicsWithoutRow()} rather than fitted.
 */
public class TopicModel {
    private final Map<String, String[]> genomePaths;
    private final Map<String, long[]> pathCounts;
    private final long threshold;
    private final double backgroundShare;
    private final int maxIterations;
    private final double tolerance;

    /**
     * Creates the model for one database.
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
        if (threshold < 1) {
            throw new IllegalArgumentException("The threshold must be at least 1: " + threshold);
        }
        if (backgroundShare < 0 || backgroundShare >= 1) {
            throw new IllegalArgumentException("The background share must be in [0, 1): " + backgroundShare);
        }
        if (maxIterations < 1) {
            throw new IllegalArgumentException("At least one iteration is needed: " + maxIterations);
        }
        for (Map.Entry<String, long[]> e : pathCounts.entrySet()) {
            String[] path = genomePaths.get(e.getKey());
            if (path == null || path.length != e.getValue().length) {
                throw new IllegalArgumentException("Path and counts of genome " + e.getKey() + " do not match.");
            }
        }
        this.genomePaths = genomePaths;
        this.pathCounts = pathCounts;
        this.threshold = threshold;
        this.backgroundShare = backgroundShare;
        this.maxIterations = maxIterations;
        this.tolerance = tolerance;
    }

    /**
     * Fits the model to one sample.
     *
     * @param nodeCounts the sample's count c<sub>n</sub> per node tax id, in k-mers
     * @return the topics and how the fit went
     */
    public TopicModelResult fit(Map<String, Long> nodeCounts) {
        // Equation (2): the admitted nodes, in a fixed order so that the fit is reproducible.
        Map<String, Integer> admitted = new LinkedHashMap<>();
        List<String> admittedNodes = new ArrayList<>();
        for (Map.Entry<String, Long> e : nodeCounts.entrySet()) {
            if (e.getKey() != null && e.getValue() != null && e.getValue() >= threshold) {
                admittedNodes.add(e.getKey());
            }
        }
        admittedNodes.sort(null);
        for (String node : admittedNodes) {
            admitted.put(node, admitted.size());
        }

        // Equation (3): each genome to the most specific admitted node on its path.
        Map<String, List<Member>> pools = new HashMap<>();
        Map<String, String[]> topicPaths = new HashMap<>();
        List<String> genomes = new ArrayList<>(pathCounts.keySet());
        genomes.sort(null);
        for (String genome : genomes) {
            String[] path = genomePaths.get(genome);
            for (int i = 0; i < path.length; i++) {
                if (admitted.containsKey(path[i])) {
                    pools.computeIfAbsent(path[i], k -> new ArrayList<>()).add(new Member(pathCounts.get(genome), i));
                    topicPaths.putIfAbsent(path[i], Arrays.copyOfRange(path, i, path.length));
                    break;
                }
            }
        }

        // Equations (4) and (5): one smoothed row per topic, kept only at the admitted nodes, since the
        // likelihood of Equation (1) runs over those alone. Topics are ordered like the admitted nodes.
        List<String> topicNodes = new ArrayList<>();
        List<String> topicsWithoutRow = new ArrayList<>();
        List<int[]> rowNodes = new ArrayList<>();
        List<double[]> rowValues = new ArrayList<>();
        List<double[]> holds = new ArrayList<>();
        for (String node : admittedNodes) {
            List<Member> pool = pools.get(node);
            String[] path = topicPaths.get(node);
            if (pool == null) {
                topicsWithoutRow.add(node);
                continue;
            }
            long[] row = new long[path.length];
            double[] memberHolds = new double[pool.size()];
            for (int j = 0; j < pool.size(); j++) {
                Member member = pool.get(j);
                for (int k = 0; k < row.length; k++) {
                    row[k] += member.counts[member.offset + k];
                    memberHolds[j] += member.counts[member.offset + k];
                }
            }
            long total = 0;
            for (long c : row) {
                total += c;
            }
            if (total == 0) {
                topicsWithoutRow.add(node);
                continue;
            }
            int onAdmitted = 0;
            for (String n : path) {
                if (admitted.containsKey(n)) {
                    onAdmitted++;
                }
            }
            int[] nodes = new int[onAdmitted];
            double[] values = new double[onAdmitted];
            int j = 0;
            for (int k = 0; k < path.length; k++) {
                Integer index = admitted.get(path[k]);
                if (index != null) {
                    nodes[j] = index;
                    values[j] = (1 - backgroundShare) * row[k] / total + backgroundShare / path.length;
                    j++;
                }
            }
            topicNodes.add(node);
            rowNodes.add(nodes);
            rowValues.add(values);
            holds.add(memberHolds);
        }

        // The counts the fit explains: those at admitted nodes some topic's row gives mass. Lying on a
        // topic's lineage is not enough without smoothing - a row can be zero there, and the node's
        // mixture with it, which the update below would divide by.
        boolean[] covered = new boolean[admittedNodes.size()];
        for (int m = 0; m < rowNodes.size(); m++) {
            int[] nodes = rowNodes.get(m);
            double[] values = rowValues.get(m);
            for (int j = 0; j < nodes.length; j++) {
                if (values[j] > 0) {
                    covered[nodes[j]] = true;
                }
            }
        }
        long[] counts = new long[admittedNodes.size()];
        long fittedCount = 0;
        List<String> uncoveredNodes = new ArrayList<>();
        for (int n = 0; n < counts.length; n++) {
            if (covered[n]) {
                counts[n] = nodeCounts.get(admittedNodes.get(n));
                fittedCount += counts[n];
            } else {
                uncoveredNodes.add(admittedNodes.get(n));
            }
        }

        // Equation (13), started uniform since the update keeps a zero share at zero.
        int topicCount = topicNodes.size();
        double[] theta = new double[topicCount];
        Arrays.fill(theta, 1.0 / topicCount);
        double[] mixture = new double[counts.length];
        double[] next = new double[topicCount];
        int iterations = 0;
        // With nothing to fit, the shares stay uniform rather than turning into 0 / 0.
        boolean converged = topicCount <= 1 || fittedCount == 0;
        while (!converged && iterations < maxIterations) {
            iterations++;
            Arrays.fill(mixture, 0);
            for (int m = 0; m < topicCount; m++) {
                int[] nodes = rowNodes.get(m);
                double[] values = rowValues.get(m);
                for (int j = 0; j < nodes.length; j++) {
                    mixture[nodes[j]] += values[j] * theta[m];
                }
            }
            double change = 0;
            for (int m = 0; m < topicCount; m++) {
                int[] nodes = rowNodes.get(m);
                double[] values = rowValues.get(m);
                double sum = 0;
                for (int j = 0; j < nodes.length; j++) {
                    if (counts[nodes[j]] > 0) {
                        sum += counts[nodes[j]] * values[j] / mixture[nodes[j]];
                    }
                }
                next[m] = theta[m] * sum / fittedCount;
                change += Math.abs(next[m] - theta[m]);
            }
            double[] swap = theta;
            theta = next;
            next = swap;
            converged = change < tolerance;
        }

        // Equations (14) to (18): the share over what one member is expected to contribute.
        double[] expectedHold = new double[topicCount];
        double abundanceSum = 0;
        double[] abundance = new double[topicCount];
        for (int m = 0; m < topicCount; m++) {
            double sum = 0;
            double squares = 0;
            for (double h : holds.get(m)) {
                sum += h;
                squares += h * h;
            }
            expectedHold[m] = squares / sum;
            abundance[m] = theta[m] / expectedHold[m];
            abundanceSum += abundance[m];
        }

        List<Topic> topics = new ArrayList<>();
        for (int m = 0; m < topicCount; m++) {
            String node = topicNodes.get(m);
            topics.add(new Topic(node, nodeCounts.get(node), pools.get(node).size(), theta[m], expectedHold[m],
                    abundanceSum > 0 ? abundance[m] / abundanceSum : 0));
        }
        topics.sort((a, b) -> Double.compare(b.getKmerShare(), a.getKmerShare()));
        return new TopicModelResult(topics, iterations, converged, fittedCount, topicsWithoutRow, uncoveredNodes);
    }

    /**
     * A genome as a member of a topic's pool: its counts and where the topic's node sits on its path.
     */
    private static class Member {
        private final long[] counts;
        private final int offset;

        private Member(long[] counts, int offset) {
            this.counts = counts;
            this.offset = offset;
        }
    }
}
