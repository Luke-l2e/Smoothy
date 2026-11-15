package org.metagene.genestrip.smoothy;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.make.ConfigKey;
import org.metagene.genestrip.make.ConfigParamInfo;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.MDDescription;
import org.metagene.genestrip.tax.Rank;

/**
 * Configuration keys added by the smoothy extension, all named with the prefix {@code smoothy.}.
 * Constants annotated with {@link MDDescription} contribute to the generated documentation.
 */
public enum SmoothyConfigKey implements ConfigKey {
    /**
     * Whether the path counts also cover the genomes the database did not store, the second and third
     * case of Section 3.2 (p. 21): a genome of a taxon the database names but whose genome it did not
     * enter, and a genome of a taxon the database's tree does not have at all.
     */
    @MDDescription("Whether `kmergenomepathcounts` also counts the genomes the database did not store, reading every "
            + "genome of the RefSeq release that took part in the database's LCA update (`refseq.updateScope`). A "
            + "genome of a taxon the database names but did not enter is filed at that taxon, one of a taxon missing "
            + "from the database's tree at the nearest ancestor the tree has. Either joins the pool of the most "
            + "specific topic above it, which gives the topics at inner nodes rows for material the database cannot "
            + "name. `false` counts the stored genomes only, as before. The setting is part of the output files' "
            + "names, so counts of both kinds can sit side by side.")
    UNSTORED_GENOMES("unstoredGenomes", new ConfigParamInfo.BooleanConfigParamInfo(true),
            SmoothyGoalKey.KMER_GENOME_PATH_COUNTS, SmoothyGoalKey.KMER_GENOME_PATH_COUNTS_CSV,
            SmoothyGoalKey.KMER_GENOME_PATH_COUNTS_SER, SmoothyGoalKey.LOAD_KMER_GENOME_PATH_COUNTS,
            SmoothyGoalKey.FT_KMER_GENOME_PATH_COUNTS, SmoothyGoalKey.FT_KMER_GENOME_PATH_COUNTS_CSV,
            SmoothyGoalKey.FT_KMER_GENOME_PATH_COUNTS_SER, SmoothyGoalKey.LOAD_FT_KMER_GENOME_PATH_COUNTS),
    /**
     * The admission threshold T of Equation (2): the count a node needs to be a topic at all.
     */
    @MDDescription("The admission threshold *T* of the topic model: a node is a topic for a sample only if the sample "
            + "reports at least this many *k*-mers at it. Counts at nodes below the threshold are left out of the fit "
            + "rather than folded into the node above.")
    TOPIC_THRESHOLD("topicThreshold", new ConfigParamInfo.LongConfigParamInfo(1, Long.MAX_VALUE, 1),
            SmoothyGoalKey.TOPIC_MODEL, SmoothyGoalKey.FT_TOPIC_MODEL),
    /**
     * The share epsilon of Equation (5) every topic's row hands to a uniform background over its path.
     */
    @MDDescription("The share *&epsilon;* of each topic's distribution that is spread evenly over the nodes from the topic's "
            + "node up to the root, so that a topic can explain some of the count at every node of its lineage.")
    BACKGROUND_SHARE("backgroundShare", new ConfigParamInfo.DoubleConfigParamInfo(0, 1, 0.05, true),
            SmoothyGoalKey.TOPIC_MODEL, SmoothyGoalKey.FT_TOPIC_MODEL),
    /**
     * The most iterations the EM fit of Equation (13) may take.
     */
    @MDDescription("The maximum number of iterations of the EM fit of the topic model. A fit that does not converge "
            + "within them is reported, and its last estimate is used.")
    MAX_ITERATIONS("maxIterations", new ConfigParamInfo.IntConfigParamInfo(1, Integer.MAX_VALUE, 10000),
            SmoothyGoalKey.TOPIC_MODEL, SmoothyGoalKey.FT_TOPIC_MODEL),
    /**
     * The change below which the EM fit counts as converged.
     */
    @MDDescription("The EM fit of the topic model stops once the topic shares change by less than this in total "
            + "(summed absolute change) from one iteration to the next.")
    CONVERGENCE_TOLERANCE("convergenceTolerance", new ConfigParamInfo.DoubleConfigParamInfo(0, 1, 1e-10, true),
            SmoothyGoalKey.TOPIC_MODEL, SmoothyGoalKey.FT_TOPIC_MODEL),
    /**
     * The severity beta of Equation (10): how hard the consistency weight attenuates an inconsistent node.
     */
    @MDDescription("The severity *&beta;* of the *k*-mer consistency weight of the topic model. A node's count enters the "
            + "fit multiplied by max(0, 1 + *&beta;* ln min(*r*, 1)), where *r* is the node's `unique kmers / exp.` "
            + "as Genestrip's match result reports it. So a node showing fewer distinct *k*-mers than uniform "
            + "sampling would give counts less, and one below *r* = e<sup>-1/&beta;</sup> not at all. `0` turns "
            + "the weight off.")
    CONSISTENCY_SEVERITY("consistencySeverity", new ConfigParamInfo.DoubleConfigParamInfo(0, 100, 0.5),
            SmoothyGoalKey.TOPIC_MODEL, SmoothyGoalKey.FT_TOPIC_MODEL),
    /**
     * The reach threshold tau of Equation (9): the reach a node needs for the consistency weight to apply.
     */
    @MDDescription("The reach *&tau;* a node needs for the *k*-mer consistency weight to be applied to it. The reach is "
            + "the least share of the node's database *k*-mers that any genome the node speaks for holds there. "
            + "Below it, the node's *k*-mers are shared by organisms of which a sample may hold only some, so a "
            + "low consistency ratio would not mean an inconsistent match. Such a node keeps its count as it is.")
    CONSISTENCY_REACH_THRESHOLD("consistencyReachThreshold", new ConfigParamInfo.DoubleConfigParamInfo(0, 1, 0.9),
            SmoothyGoalKey.TOPIC_MODEL, SmoothyGoalKey.FT_TOPIC_MODEL),
    /**
     * The member floor gamma of Equation (8): the pool weight a genome needs to enter the reach.
     */
    @MDDescription("The share of a node's pool *k*-mers a genome needs to count as a member when the reach is taken. "
            + "Smaller genomes, such as fragments of a draft assembly beside a finished genome, are left out, so "
            + "that they do not decline the weight for the whole node.")
    CONSISTENCY_MEMBER_FLOOR("consistencyMemberFloor", new ConfigParamInfo.DoubleConfigParamInfo(0, 1, 0.01),
            SmoothyGoalKey.TOPIC_MODEL, SmoothyGoalKey.FT_TOPIC_MODEL),
    /**
     * The rank gate rho of Equation (9): the consistency weight is considered at this rank and below only.
     */
    @MDDescription("The rank *&rho;* at and below which the *k*-mer consistency weight is considered at all; a node "
            + "without a rank of its own counts as below *&rho;* if one of its ancestors is at or below it. Higher "
            + "nodes keep their counts as they are.")
    CONSISTENCY_RANK("consistencyRank", new GSConfigKey.RankConfigParamInfo(Rank.SPECIES),
            SmoothyGoalKey.TOPIC_MODEL, SmoothyGoalKey.FT_TOPIC_MODEL);

    /**
     * The prefix of every smoothy configuration key.
     */
    public static final String SMOOTHY = "smoothy";

    private final String name;
    private final ConfigParamInfo<?> param;
    private final SmoothyGoalKey[] forGoals;

    SmoothyConfigKey(String name, ConfigParamInfo<?> param, SmoothyGoalKey... forGoals) {
        this.name = SMOOTHY + "." + name;
        this.param = param;
        this.forGoals = forGoals;
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public ConfigParamInfo<?> getInfo() {
        return param;
    }

    @Override
    public boolean isForGoal(GoalKey forGoal) {
        if (forGoal == null) {
            return true;
        }
        for (GoalKey id : forGoals) {
            if (forGoal.equals(id)) {
                return true;
            }
        }
        return false;
    }

    @Override
    public String toString() {
        return getName();
    }
}
