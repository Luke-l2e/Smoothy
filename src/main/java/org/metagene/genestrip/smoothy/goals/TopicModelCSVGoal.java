package org.metagene.genestrip.smoothy.goals;

import org.metagene.genestrip.GSProject.GSFileType;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.smoothy.SmoothyProject;
import org.metagene.genestrip.smoothy.model.Topic;
import org.metagene.genestrip.smoothy.model.TopicModelResult;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Writes the topics {@link TopicModelGoal} fitted to one CSV file per sample, named and placed like the
 * files of Genestrip's {@code match} goal: in the project's {@code csv} folder, with the sample's key in
 * the name.
 * <p>
 * One row per topic, by descending k-mer share. The columns are:
 * <ul>
 * <li>{@code taxid}, {@code name} and {@code rank} - the topic's node; at a genome node the name is the
 * genome's key, e.g. {@code NC_001477};</li>
 * <li>{@code kmers} - the count c<sub>m</sub> the sample reports at the node;</li>
 * <li>{@code genomes} - how many genomes the topic speaks for: one at a genome node, more at a taxon,
 * whose topic is then a share of the clade and not of one organism;</li>
 * <li>{@code kmer share} - &theta;<sub>m</sub>, the share of the sample's k-mers the topic accounts for;</li>
 * <li>{@code expected hold} - E<sub>m</sub>, the k-mers one organism of the topic is expected to
 * contribute;</li>
 * <li>{@code abundance} - the topic's relative abundance, its share of organisms rather than of
 * sequence.</li>
 * </ul>
 *
 * @param <P> the project type
 */
public class TopicModelCSVGoal<P extends SmoothyProject> extends FileListGoal<P> {
    private static final DecimalFormat DF = new DecimalFormat("0.00000000", new DecimalFormatSymbols(Locale.US));

    private final ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal;
    private final ObjectGoal<Map<String, TopicModelResult>, P> topicModelGoal;
    private final ObjectGoal<Database, P> dbGoal;
    private final Map<File, String> fileToKeyMap;

    /**
     * Creates the goal.
     *
     * @param project        the project this goal belongs to
     * @param key            the key identifying this goal
     * @param fastqMapGoal   the goal supplying the samples, whose keys name the files
     * @param topicModelGoal the goal fitting the topic model per sample
     * @param dbGoal         the goal supplying the database, whose tree names the topics' nodes
     * @param deps           any further goals this goal depends on
     */
    @SafeVarargs
    public TopicModelCSVGoal(P project, GoalKey key, ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal,
                             ObjectGoal<Map<String, TopicModelResult>, P> topicModelGoal,
                             ObjectGoal<Database, P> dbGoal, Goal<P>... deps) {
        super(project, key, (List<File>) null, Goal.append(deps, fastqMapGoal, topicModelGoal, dbGoal));
        this.fastqMapGoal = fastqMapGoal;
        this.topicModelGoal = topicModelGoal;
        this.dbGoal = dbGoal;
        fileToKeyMap = new HashMap<>();
    }

    @Override
    protected void provideFiles() {
        for (String sample : fastqMapGoal.get().keySet()) {
            File file = getProject().getOutputFile(getKey().getName(), sample, null, GSFileType.CSV, false);
            addFile(file);
            fileToKeyMap.put(file, sample);
        }
    }

    @Override
    protected void makeFile(File file) throws IOException {
        TopicModelResult result = topicModelGoal.get().get(fileToKeyMap.get(file));
        SmallTaxTree tree = dbGoal.get().getTaxTree();
        try (PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file), false,
                StandardCharsets.UTF_8.name())) {
            out.println("taxid;name;rank;kmers;genomes;kmer share;expected hold;abundance;");
            for (Topic topic : result.getTopics()) {
                SmallTaxIdNode node = tree.getNodeByTaxId(topic.getNode());
                Rank rank = node == null ? null : node.getRank();
                out.print(topic.getNode());
                out.print(';');
                out.print(node == null ? "" : node.getName());
                out.print(';');
                out.print(rank == null ? "" : rank.getName());
                out.print(';');
                out.print(topic.getCount());
                out.print(';');
                out.print(topic.getGenomes());
                out.print(';');
                out.print(DF.format(topic.getKmerShare()));
                out.print(';');
                out.print(DF.format(topic.getExpectedHold()));
                out.print(';');
                out.print(DF.format(topic.getAbundance()));
                out.println(';');
            }
        }
    }
}
