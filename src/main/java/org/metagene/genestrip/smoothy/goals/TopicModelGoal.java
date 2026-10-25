package org.metagene.genestrip.smoothy.goals;

import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.CountsPerTaxid;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.smoothy.SmoothyConfigKey;
import org.metagene.genestrip.smoothy.SmoothyProject;
import org.metagene.genestrip.smoothy.model.TopicModel;
import org.metagene.genestrip.smoothy.model.TopicModelResult;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import java.util.*;

/**
 * Fits the {@link TopicModel} to every sample of the project, giving per sample the share of its k-mers
 * and the relative abundance of each topic (Sections 3.1 to 3.6 of the smoothy paper).
 * <p>
 * The two inputs are the per-genome path counts, supplied by {@link LoadKMerGenomePathCountsGoal} so
 * that they are counted once per database and not once per run, and the per-node k-mer counts
 * c<sub>n</sub> the match result reports for each sample - {@link CountsPerTaxid#getKMers()}, the k-mers
 * and not the reads, for the reasons Section 3.1 gives. Both have to come from the same database: the
 * paths are read off its tree, and its nodes are what the counts of either side refer to. The goal is
 * therefore registered once for Genestrip's database and once for the FT database, each with the match
 * result and the path counts of its own database.
 * <p>
 * The result maps each sample key, as the match result has it, to the fitted model.
 *
 * @param <P> the project type
 */
public class TopicModelGoal<P extends SmoothyProject> extends ObjectGoal<Map<String, TopicModelResult>, P> {
    private final ObjectGoal<Map<String, long[]>, P> pathCountsGoal;
    private final ObjectGoal<Map<String, MatchingResult>, P> matchResGoal;
    private final ObjectGoal<Database, P> dbGoal;

    /**
     * Creates the goal.
     *
     * @param project        the project this goal belongs to
     * @param key            the key identifying this goal
     * @param pathCountsGoal the goal supplying the per-genome path counts of the database
     * @param matchResGoal   the goal supplying the match result per sample, against the same database
     * @param dbGoal         the goal supplying that database, whose tree gives the genomes' paths
     * @param deps           any further goals this goal depends on
     */
    @SafeVarargs
    public TopicModelGoal(P project, GoalKey key, ObjectGoal<Map<String, long[]>, P> pathCountsGoal,
                          ObjectGoal<Map<String, MatchingResult>, P> matchResGoal, ObjectGoal<Database, P> dbGoal,
                          Goal<P>... deps) {
        super(project, key, Goal.append(deps, pathCountsGoal, matchResGoal, dbGoal));
        this.pathCountsGoal = pathCountsGoal;
        this.matchResGoal = matchResGoal;
        this.dbGoal = dbGoal;
    }

    @Override
    protected void doMakeThis() {
        Map<String, long[]> pathCounts = pathCountsGoal.get();
        TopicModel model = new TopicModel(genomePaths(dbGoal.get().getTaxTree(), pathCounts), pathCounts,
                longConfigValue(SmoothyConfigKey.TOPIC_THRESHOLD), doubleConfigValue(SmoothyConfigKey.BACKGROUND_SHARE),
                intConfigValue(SmoothyConfigKey.MAX_ITERATIONS), doubleConfigValue(SmoothyConfigKey.CONVERGENCE_TOLERANCE));

        Map<String, TopicModelResult> results = new LinkedHashMap<>(); // In the order of the samples.
        for (Map.Entry<String, MatchingResult> e : matchResGoal.get().entrySet()) {
            TopicModelResult result = model.fit(nodeCounts(e.getValue()));
            report(e.getKey(), result);
            results.put(e.getKey(), result);
        }
        set(results);
    }

    /**
     * Reads each genome's path to the root off the database's tree, as tax ids from the genome's own
     * node up.
     *
     * @param tree       the tree of the database the counts were measured on
     * @param pathCounts the per-genome path counts, whose keys are the genomes' tax ids
     * @return per genome, the tax ids along its path
     * @throws IllegalStateException if a genome is not in the tree or its path does not match its counts,
     *                               i.e. counts and tree belong to different databases
     */
    protected static Map<String, String[]> genomePaths(SmallTaxTree tree, Map<String, long[]> pathCounts) {
        Map<String, String[]> paths = new HashMap<>();
        for (Map.Entry<String, long[]> e : pathCounts.entrySet()) {
            SmallTaxIdNode genome = tree.getNodeByTaxId(e.getKey());
            if (genome == null) {
                throw new IllegalStateException("Genome " + e.getKey() + " of the path counts is not in the database's tree.");
            }
            List<String> path = new ArrayList<>();
            for (SmallTaxIdNode n = genome; n != null; n = n.getParent()) {
                path.add(n.getTaxId());
            }
            if (path.size() != e.getValue().length) {
                throw new IllegalStateException("The path of genome " + e.getKey() + " in the database's tree does not"
                        + " match its counts: " + path.size() + " nodes for " + e.getValue().length + " counts.");
            }
            paths.put(e.getKey(), path.toArray(new String[0]));
        }
        return paths;
    }

    /**
     * Returns the sample's k-mer count per node, c<sub>n</sub> of the paper.
     *
     * @param matchingResult the sample's match result
     * @return the k-mers reported per node tax id
     */
    protected static Map<String, Long> nodeCounts(MatchingResult matchingResult) {
        Map<String, Long> counts = new HashMap<>();
        for (Map.Entry<String, CountsPerTaxid> e : matchingResult.getTaxid2Stats().entrySet()) {
            // The entry without a tax id holds the totals over all nodes, not a node of its own.
            if (e.getKey() != null) {
                counts.put(e.getKey(), e.getValue().getKMers());
            }
        }
        return counts;
    }

    private void report(String sample, TopicModelResult result) {
        if (getLogger().isInfoEnabled()) {
            getLogger().info("Sample " + sample + ": " + result.getTopics().size() + " topics fitted to "
                    + result.getFittedCount() + " k-mers in " + result.getIterations() + " iterations.");
        }
        if (!result.isConverged() && getLogger().isWarnEnabled()) {
            getLogger().warn("Sample " + sample + ": the fit did not converge within "
                    + intConfigValue(SmoothyConfigKey.MAX_ITERATIONS) + " iterations, which near-identical topics"
                    + " cause. Its last estimate is used.");
        }
        if (!result.getUncoveredNodes().isEmpty() && getLogger().isWarnEnabled()) {
            getLogger().warn("Sample " + sample + ": no topic reaches the nodes " + result.getUncoveredNodes()
                    + ", so their counts are left out. Do the match result and the path counts belong to the same"
                    + " database?");
        }
        if (!result.getTopicsWithoutRow().isEmpty() && getLogger().isInfoEnabled()) {
            getLogger().info("Sample " + sample + ": " + result.getTopicsWithoutRow().size() + " admitted nodes are no"
                    + " topics, every genome beneath them having a topic of its own.");
        }
    }
}
