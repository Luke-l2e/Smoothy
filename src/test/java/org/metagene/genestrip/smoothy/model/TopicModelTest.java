package org.metagene.genestrip.smoothy.model;

import org.junit.Test;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.*;

/**
 * Tests of {@link TopicModel} on the tree of Figure 1 of the smoothy paper, built by hand:
 *
 * <pre>
 *          r
 *        /   \
 *      n1     n2
 *     /  \     \
 *   g1    g2    g3
 * </pre>
 *
 * g1 and g2 share k-mers at n1, all three at r. The rows of the three genomes are
 * g1: (60, 30, 10), g2: (50, 40, 10) and g3: (180, 0, 20) on their paths, so g3 holds nothing at n2 and
 * is twice the size of the others.
 */
public class TopicModelTest {
    private static final double EXACT = 1e-6;

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

    /**
     * A sample drawn exactly from the mixture 0.5 g1 + 0.3 g2 + 0.2 g3, i.e. counts
     * 1000 * (0.5 * row(g1) + 0.3 * row(g2) + 0.2 * row(g3)) with the rows normalised.
     */
    private static Map<String, Long> mixtureSample() {
        Map<String, Long> c = new HashMap<>();
        c.put("g1", 300L);
        c.put("g2", 150L);
        c.put("g3", 180L);
        c.put("n1", 270L);
        c.put("r", 100L);
        return c;
    }

    private static TopicModel model(long threshold, double backgroundShare) {
        return new TopicModel(paths(), pathCounts(), threshold, backgroundShare, 100000, 1e-13);
    }

    private static Topic topic(TopicModelResult result, String node) {
        for (Topic t : result.getTopics()) {
            if (t.getNode().equals(node)) {
                return t;
            }
        }
        return null;
    }

    /** Without smoothing, counts drawn exactly from a mixture give back exactly that mixture. */
    @Test
    public void exactMixtureIsRecovered() {
        TopicModelResult result = model(1, 0).fit(mixtureSample());

        assertTrue(result.isConverged());
        assertEquals(0.5, topic(result, "g1").getKmerShare(), EXACT);
        assertEquals(0.3, topic(result, "g2").getKmerShare(), EXACT);
        assertEquals(0.2, topic(result, "g3").getKmerShare(), EXACT);
        // The topics come by descending share.
        assertEquals("g1", result.getTopics().get(0).getNode());
    }

    /**
     * n1 and r pass the threshold but every genome beneath them has a topic of its own, so no genome
     * gives them a row: they are reported, not fitted, and their counts are still explained by the
     * genomes' rows.
     */
    @Test
    public void taxaWhoseGenomesAllHaveTopicsGetNoRow() {
        TopicModelResult result = model(1, 0).fit(mixtureSample());

        assertEquals(3, result.getTopics().size());
        assertTrue(result.getTopicsWithoutRow().contains("n1"));
        assertTrue(result.getTopicsWithoutRow().contains("r"));
        assertTrue(result.getUncoveredNodes().isEmpty());
        assertEquals(1000, result.getFittedCount());
    }

    /**
     * Abundance divides each share by what one member contributes - here a genome's whole size - so
     * g3, twice as large, has half the abundance its k-mer share suggests.
     */
    @Test
    public void abundanceDividesByGenomeSize() {
        TopicModelResult result = model(1, 0).fit(mixtureSample());

        assertEquals(100, topic(result, "g1").getExpectedHold(), EXACT);
        assertEquals(200, topic(result, "g3").getExpectedHold(), EXACT);
        double norm = 0.5 / 100 + 0.3 / 100 + 0.2 / 200;
        assertEquals(0.5 / 100 / norm, topic(result, "g1").getAbundance(), EXACT);
        assertEquals(0.2 / 200 / norm, topic(result, "g3").getAbundance(), EXACT);
    }

    /**
     * A genome without a count at its own node is no topic; it falls to the nearest topic above it, and
     * only what it holds from that node upwards enters the row there (Figure 3 of the paper).
     */
    @Test
    public void genomeWithoutTopicFallsToTheNodeAbove() {
        Map<String, Long> c = new HashMap<>();
        c.put("g1", 60L);
        c.put("n1", 30L);
        c.put("r", 10L);
        TopicModelResult result = model(1, 0.05).fit(c);

        // g2 has no count of its own and so belongs to n1, g3 to r.
        assertEquals(1, topic(result, "n1").getGenomes());
        assertEquals(1, topic(result, "r").getGenomes());
        assertEquals(1, topic(result, "g1").getGenomes());
        // What one member of n1 holds is g2 from n1 upwards: 40 + 10, not its 50 at g2.
        assertEquals(50, topic(result, "n1").getExpectedHold(), EXACT);
        assertEquals(20, topic(result, "r").getExpectedHold(), EXACT);
        // The sample is g1 alone, so g1 carries nearly all of it.
        assertTrue(topic(result, "g1").getKmerShare() > 0.9);
    }

    /** Counts below the threshold are left out of the fit rather than folded into the node above. */
    @Test
    public void countsBelowTheThresholdAreLeftOut() {
        TopicModelResult result = model(200, 0.05).fit(mixtureSample());

        // Only g1 (300) and n1 (270) pass; g2 falls to n1, g3 has no admitted node on its path.
        assertEquals(2, result.getTopics().size());
        assertEquals(1, topic(result, "n1").getGenomes());
        assertNull(topic(result, "g3"));
        assertEquals(570, result.getFittedCount());
    }

    /** A count at a node no genome's path reaches cannot be explained by any topic. */
    @Test
    public void countOffEveryPathIsReported() {
        Map<String, Long> c = mixtureSample();
        c.put("x", 5L);
        TopicModelResult result = model(1, 0.05).fit(c);

        assertEquals(List.of("x"), result.getUncoveredNodes());
        assertEquals(1000, result.getFittedCount());
    }

    /**
     * Without smoothing, a node on a topic's lineage where its row is zero cannot be explained by it:
     * g3 holds nothing at n2. The count there is reported and left out instead of dividing by zero.
     */
    @Test
    public void zeroRowWithoutSmoothingLeavesTheNodeOut() {
        Map<String, Long> c = new HashMap<>();
        c.put("g3", 90L);
        c.put("g1", 60L);
        c.put("n2", 5L);
        TopicModelResult result = model(1, 0).fit(c);

        assertEquals(List.of("n2"), result.getUncoveredNodes());
        assertEquals(150, result.getFittedCount());
        for (Topic t : result.getTopics()) {
            assertFalse(t.getNode(), Double.isNaN(t.getKmerShare()));
        }
    }

    /** With smoothing the shares are still a distribution over the topics. */
    @Test
    public void sharesSumToOne() {
        TopicModelResult result = model(1, 0.05).fit(mixtureSample());

        double kmerShares = 0;
        double abundances = 0;
        for (Topic t : result.getTopics()) {
            assertTrue(t.getKmerShare() >= 0);
            kmerShares += t.getKmerShare();
            abundances += t.getAbundance();
        }
        assertEquals(1, kmerShares, EXACT);
        assertEquals(1, abundances, EXACT);
    }

    /** A sample with a single topic needs no iteration at all. */
    @Test
    public void singleTopicTakesEverything() {
        Map<String, Long> c = new HashMap<>();
        c.put("g3", 90L);
        TopicModelResult result = model(1, 0.05).fit(c);

        assertEquals(1, result.getTopics().size());
        assertEquals(1, topic(result, "g3").getKmerShare(), EXACT);
        assertEquals(0, result.getIterations());
    }
}
