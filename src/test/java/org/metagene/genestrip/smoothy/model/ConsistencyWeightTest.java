package org.metagene.genestrip.smoothy.model;

import org.junit.Test;

import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.*;

/**
 * Tests of the k-mer consistency weight of Section 3.3 of the smoothy paper, Equations (7) to (10), alone
 * and inside {@link TopicModel#fit(Map, Map)}. The tree is that of {@link TopicModelTest}:
 *
 * <pre>
 *          r
 *        /   \
 *      n1     n2
 *     /  \     \
 *   g1    g2    g3
 * </pre>
 *
 * with the rows g1: (60, 30, 10), g2: (50, 40, 10) and g3: (180, 0, 20) on their paths.
 */
public class ConsistencyWeightTest {
    private static final double EXACT = 1e-9;
    private static final ConsistencyWeight PAPER = new ConsistencyWeight(0.5, 0.9, 0.01); // Table 5, γ chosen

    private static Map<String, String[]> paths() {
        Map<String, String[]> paths = new HashMap<>();
        paths.put("g1", new String[]{"g1", "n1", "r"});
        paths.put("g2", new String[]{"g2", "n1", "r"});
        paths.put("g3", new String[]{"g3", "n2", "r"});
        return paths;
    }

    private static Map<String, long[]> pathCounts() {
        Map<String, long[]> counts = new HashMap<>();
        counts.put("g1", new long[]{60, 30, 10});
        counts.put("g2", new long[]{50, 40, 10});
        counts.put("g3", new long[]{180, 0, 20});
        return counts;
    }

    /** The sample 0.5 g1 + 0.3 g2 + 0.2 g3 of {@link TopicModelTest}, 1000 k-mers. */
    private static Map<String, Long> mixtureSample() {
        Map<String, Long> c = new HashMap<>();
        c.put("g1", 300L);
        c.put("g2", 150L);
        c.put("g3", 180L);
        c.put("n1", 270L);
        c.put("r", 100L);
        return c;
    }

    private static TopicModel model(ConsistencyWeight weight) {
        return new TopicModel(paths(), pathCounts(), 1, 0, 100000, 1e-13, weight);
    }

    private static Topic topic(TopicModelResult result, String node) {
        for (Topic t : result.getTopics()) {
            if (t.getNode().equals(node)) {
                return t;
            }
        }
        return null;
    }

    /** Equation (7): E = u (1 - (1 - 1/u)^c), also for a large u where the naive form loses digits. */
    @Test
    public void expectedUniqueKMersFollowsEquation7() {
        assertEquals(60 * (1 - Math.pow(1 - 1d / 60, 300)), ConsistencyWeight.expectedUniqueKMers(60, 300), EXACT);
        assertEquals(1, ConsistencyWeight.expectedUniqueKMers(1, 5), EXACT);
        // Few draws from a huge pool hardly ever repeat.
        assertEquals(1000, ConsistencyWeight.expectedUniqueKMers(1_000_000_000_000L, 1000), 1e-3);
    }

    /** Equation (10) and the requirements listed after it. */
    @Test
    public void weightFollowsEquation10() {
        assertEquals(1, PAPER.weight(1), EXACT); // w(1) = 1
        assertEquals(1, PAPER.weight(1.3), EXACT); // no information above one
        assertEquals(1, PAPER.weight(Double.NaN), EXACT);
        // The calibration point of the paper: r = 0.8631 keeps 0.93 under β = 1/2.
        assertEquals(1 + 0.5 * Math.log(0.8631), PAPER.weight(0.8631), EXACT);
        assertEquals(0.93, PAPER.weight(0.8631), 0.005);
        // Zero from r = e^(-1/β) down.
        assertEquals(0, PAPER.weight(Math.exp(-2)), EXACT);
        assertEquals(0, PAPER.weight(0.01), EXACT);
        assertFalse(new ConsistencyWeight(0, 0.9, 0.01).isActive());
    }

    /** Without consistency data, or with the weight off, the fit is the unweighted one. */
    @Test
    public void inactiveWeightChangesNothing() {
        Map<String, NodeConsistency> inconsistent = new HashMap<>();
        inconsistent.put("g1", new NodeConsistency(10, 60, true));
        TopicModelResult plain = model(PAPER).fit(mixtureSample());
        TopicModelResult off = model(ConsistencyWeight.NONE).fit(mixtureSample(), inconsistent);

        assertEquals(0.5, topic(plain, "g1").getKmerShare(), 1e-6);
        assertEquals(0.5, topic(off, "g1").getKmerShare(), 1e-6);
        assertEquals(1, topic(off, "g1").getWeight(), EXACT);
        assertTrue(off.getNodeWeights().isEmpty());
        assertEquals(off.getFittedCount(), off.getWeightedCount(), EXACT);
    }

    /**
     * A genome node reaches all it stores (f = 1), so a ratio below one attenuates its count: g1 shows
     * half the distinct k-mers expected, its count enters with weight 1 + ln(r) / 2 and its share falls.
     */
    @Test
    public void inconsistentGenomeNodeIsAttenuated() {
        long dbKMers = 60; // = κ_g1(g1): g1 holds every k-mer stored at its node once
        long count = 300;
        long unique = Math.round(ConsistencyWeight.expectedUniqueKMers(dbKMers, count) / 2);
        double ratio = ConsistencyWeight.ratio(unique, dbKMers, count);
        Map<String, NodeConsistency> consistency = new HashMap<>();
        consistency.put("g1", new NodeConsistency(unique, dbKMers, true));
        // g2 is consistent, g3 outside the rank gate although inconsistent.
        consistency.put("g2", new NodeConsistency(Math.round(ConsistencyWeight.expectedUniqueKMers(50, 150)), 50, true));
        consistency.put("g3", new NodeConsistency(1, 180, false));

        TopicModelResult result = model(PAPER).fit(mixtureSample(), consistency);

        double weight = 1 + 0.5 * Math.log(ratio);
        assertEquals(weight, topic(result, "g1").getWeight(), EXACT);
        assertEquals(1, topic(result, "g2").getWeight(), 1e-3);
        assertEquals(1, topic(result, "g3").getWeight(), EXACT);
        assertTrue(topic(result, "g1").getKmerShare() < 0.5);
        assertEquals(1000, result.getFittedCount());
        assertEquals(1000 - (1 - weight) * count, result.getWeightedCount(), 1e-3);

        assertEquals(2, result.getNodeWeights().size()); // g1 and g2; g3 is outside the gate
        NodeWeight g1 = result.getNodeWeights().get(0);
        assertEquals("g1", g1.getNode());
        assertTrue(g1.isApplied());
        assertEquals(1, g1.getReach(), EXACT);
        assertEquals(ratio, g1.getRatio(), EXACT);
    }

    /**
     * Where the database stores more at a node than any member holds there, the reach falls below τ and
     * the weight is declined: the node keeps its count, and the decision is still reported.
     */
    @Test
    public void lowReachDeclinesTheWeight() {
        Map<String, NodeConsistency> consistency = new HashMap<>();
        consistency.put("g1", new NodeConsistency(10, 100, true)); // f = 60 / 100 < 0.9

        TopicModelResult result = model(PAPER).fit(mixtureSample(), consistency);

        assertEquals(1, topic(result, "g1").getWeight(), EXACT);
        assertEquals(0.5, topic(result, "g1").getKmerShare(), 1e-6);
        NodeWeight g1 = result.getNodeWeights().get(0);
        assertFalse(g1.isApplied());
        assertEquals(0.6, g1.getReach(), EXACT);
    }

    /**
     * At n1, the topic of g1 and g2 once their own nodes are not admitted, the reach is the least share
     * of n1's stored k-mers a member holds: 30 / 40 for g1, 40 / 40 for g2. The minimum declines the
     * weight; a member floor γ above g1's pool weight q = 40 / 90 leaves g2 alone, and the weight applies.
     */
    @Test
    public void memberFloorBoundsTheMinimum() {
        Map<String, Long> sample = new HashMap<>();
        sample.put("n1", 400L);
        sample.put("r", 100L);
        Map<String, NodeConsistency> consistency = new HashMap<>();
        consistency.put("n1", new NodeConsistency(20, 40, true));

        TopicModelResult all = model(new ConsistencyWeight(0.5, 0.9, 0.01)).fit(sample, consistency);
        NodeWeight declined = all.getNodeWeights().get(0);
        assertEquals(0.75, declined.getReach(), EXACT);
        assertFalse(declined.isApplied());

        TopicModelResult floored = model(new ConsistencyWeight(0.5, 0.9, 0.5)).fit(sample, consistency);
        NodeWeight applied = floored.getNodeWeights().get(0);
        assertEquals(1, applied.getReach(), EXACT);
        assertTrue(applied.isApplied());
        assertEquals(1 + 0.5 * Math.log(ConsistencyWeight.ratio(20, 40, 400)), applied.getWeight(), EXACT);
        assertEquals(applied.getWeight(), topic(floored, "n1").getWeight(), EXACT);
    }
}
