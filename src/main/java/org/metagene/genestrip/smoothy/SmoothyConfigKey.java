package org.metagene.genestrip.smoothy;

import org.metagene.genestrip.make.ConfigKey;
import org.metagene.genestrip.make.ConfigParamInfo;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.MDDescription;

/**
 * Configuration keys added by the smoothy extension, all named with the prefix {@code smoothy.}.
 * Constants annotated with {@link MDDescription} contribute to the generated documentation.
 */
public enum SmoothyConfigKey implements ConfigKey {
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
