/*
 *
 * “Commons Clause” License Condition v1.0
 *
 * The Software is provided to you by the Licensor under the License,
 * as defined below, subject to the following condition.
 *
 * Without limiting other conditions in the License, the grant of rights under the License
 * will not include, and the License does not grant to you, the right to Sell the Software.
 *
 * For purposes of the foregoing, “Sell” means practicing any or all of the rights granted
 * to you under the License to provide to third parties, for a fee or other consideration
 * (including without limitation fees for hosting or consulting/ support services related to
 * the Software), a product or service whose value derives, entirely or substantially, from the
 * functionality of the Software. Any license notice or attribution required by the License
 * must also include this Commons Clause License Condition notice.
 *
 * Software: genestrip
 *
 * License: Apache 2.0
 *
 * Licensor: Daniel Pfeifer (daniel.pfeifer@progotec.de)
 *
 */
package org.metagene.genestrip.smoothy.goals;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.goals.refseq.FastaReaderGoal;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.*;
import org.metagene.genestrip.smoothy.SmoothyProject;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.tax.TaxNodeSelection;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Goal computing, for every genome, how many of its k-mer occurrences the database holds at each node
 * from the node the genome sits at up to the root: &kappa;<sub>h</sub>(n) of the paper.
 * <p>
 * It re-reads the reference genomes and, for every k-mer position (subject only to DUST filtering and
 * {@code kMerSampling}, exactly as {@link org.metagene.genestrip.goals.refseq.KMerRankStatsGoal} does),
 * looks the k-mer up directly in the already-built {@link Database}'s {@link KMerStore}. No Bloom
 * filter is involved and nothing is deduplicated: a k-mer occurring twice in a genome is counted twice,
 * which is what distinguishes this walk from the one that fills the database.
 * <p>
 * Section 3.2 of the paper (p. 21) tells three kinds of genome apart, and the goal counts all three:
 * <ol>
 * <li>A genome the database stored sits at its own {@code GENOME} node, which requires the database to
 * have been built with {@code genomeNodes=true}. Its entry is keyed by that node's tax id.</li>
 * <li>A genome of a taxon the database's tree does not have at all - material outside the database -
 * sits at the nearest ancestor the tree does have.</li>
 * <li>A genome of a taxon the tree has, but which the build did not enter - Genestrip enters at most
 * {@code maxGenomesPerTaxid} genomes per taxon - sits at that taxon, or at its {@code DATA} child where
 * the build gave it one, since that is where the update put the k-mers it shares with the stored
 * genomes of the taxon.</li>
 * </ol>
 * The latter two are keyed by {@link #unstoredGenomeKey(String, String)}: the tax id of the node they
 * sit at and the genome key of their accession, the same key the build would have named their
 * {@code GENOME} node with. Which genomes those are is decided as the database's LCA update decided
 * it, by {@code refseq.updateScope}, so the goal reads the RefSeq release whether or not the database
 * was filled from it. A genome none of whose k-mers is in the database gets no entry, since it adds
 * nothing to any row; a stored genome always has one.
 * <p>
 * Each value is an array indexed by position along the path from the genome's node to the root: index
 * {@code 0} is that node, and each following index is one ancestor further up. A k-mer's occurrence is
 * tallied at the index of the node the database stores it at, which - by construction of the
 * lowest-common-ancestor update the database was built with - lies on that path for every genome the
 * update took part in. The paper (p. 22) is explicit that this holds only as far as the update scope
 * reaches.
 *
 * @param <P> the project type
 */
public class KMerGenomePathCountsGoal<P extends SmoothyProject> extends FastaReaderGoal<Map<String, long[]>, P> {
    /**
     * Separates the node tax id from the genome key in the key of a genome the database did not store.
     * Tax ids are digits only, so it cannot occur in the key of a stored genome.
     */
    public static final char UNSTORED_SEPARATOR = ':';

    private final ObjectGoal<AccessionMap, P> accessionMapGoal;
    private final ObjectGoal<Database, P> dbGoal;
    private final GSConfigKey.UpdateScope updateScope;

    // Populated for the duration of doMakeThis() and read by the (multi-threaded) fasta readers.
    private KMerStore<String> store;
    private SmallTaxTree taxTree;
    private Map<String, long[]> map;
    private Set<TaxIdNode> selectedTaxNodes;
    private Set<TaxIdNode> excludedTaxNodes;
    private boolean countUnstored;
    // K-mer occurrences of unstored genomes the database holds off their path, see handleStore().
    private final AtomicLong offPath = new AtomicLong();

    /**
     * Creates the goal, wiring the accession map and the loaded database it looks k-mers up in, in
     * addition to the standard FASTA-reader dependencies.
     *
     * @param project          the project this goal belongs to
     * @param key              the key identifying this goal
     * @param bundle           the execution context supplying the worker threads
     * @param categoriesGoal   the goal supplying the RefSeq categories to read
     * @param taxNodesGoal     the goal supplying the required taxonomy nodes
     * @param fnaFilesGoal     the goal supplying the downloaded RefSeq {@code .fna} files
     * @param additionalGoal   the goal supplying additional FASTA files mapped to their tax node
     * @param accessionMapGoal the goal supplying the accession-to-taxid map
     * @param dbGoal           the goal supplying the loaded database to look k-mers up in
     * @param taxTreeGoal      the goal supplying the taxonomy tree; kept as a direct dependency so that
     *                         the (aggressively cleaned) tree is not freed and rebuilt as a second
     *                         instance while this goal runs, same reason as in {@code KMerRankStatsGoal}
     * @param deps             any further goals this goal depends on
     */
    @SafeVarargs
    public KMerGenomePathCountsGoal(P project, GoalKey key, ExecutionContext bundle,
                                    ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                                    ObjectGoal<TaxNodeSelection, P> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                                    ObjectGoal<Map<File, TaxIdNode>, P> additionalGoal,
                                    ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<Database, P> dbGoal,
                                    ObjectGoal<TaxTree, P> taxTreeGoal, Goal<P>... deps) {
        // The unstored genomes are the release's genomes the update read and the fill did not enter, so
        // counting them needs the release even where the database was not filled from it - as the
        // update itself does, see DBGoal.
        super(project, key, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal,
                true, Goal.append(deps, accessionMapGoal, dbGoal, taxTreeGoal));
        this.accessionMapGoal = accessionMapGoal;
        this.dbGoal = dbGoal;
        updateScope = (GSConfigKey.UpdateScope) project.configValue(GSConfigKey.UPDATE_SCOPE);
    }

    /**
     * Returns the key of a genome the database did not store.
     *
     * @param nodeTaxId the tax id of the node the genome sits at
     * @param genomeKey the genome key of its accession, e.g. {@code NZ_ABCD} for a WGS assembly
     * @return the key, which names both
     */
    public static String unstoredGenomeKey(String nodeTaxId, String genomeKey) {
        return nodeTaxId + UNSTORED_SEPARATOR + genomeKey;
    }

    /**
     * Whether a key of the goal's result names a genome the database did not store.
     *
     * @param key a key of the result
     * @return whether it was made by {@link #unstoredGenomeKey(String, String)}
     */
    public static boolean isUnstored(String key) {
        return key.indexOf(UNSTORED_SEPARATOR) >= 0;
    }

    /**
     * Returns the tax id of the node a genome sits at: its own {@code GENOME} node if it was stored,
     * else the node {@link #unstoredGenomeKey(String, String)} names.
     *
     * @param key a key of the result
     * @return the tax id of the node index {@code 0} of the genome's counts refers to
     */
    public static String nodeOf(String key) {
        int i = key.indexOf(UNSTORED_SEPARATOR);
        return i < 0 ? key : key.substring(0, i);
    }

    /**
     * Returns the genome key of a genome the database did not store.
     *
     * @param key a key of the result
     * @return the genome key, or {@code null} for a stored genome, whose node's name is its genome key
     */
    public static String genomeKeyOf(String key) {
        int i = key.indexOf(UNSTORED_SEPARATOR);
        return i < 0 ? null : key.substring(i + 1);
    }

    @Override
    protected void doMakeThis() {
        try {
            // Concurrent, since the entries of unstored genomes are only created when a reader first
            // finds one of their k-mers in the database. The entries of stored genomes are created here
            // beforehand, so that a stored genome always has one, even if it held nothing.
            map = new ConcurrentHashMap<>();
            Database database = dbGoal.get();
            store = database.getKmerStore();
            taxTree = database.getTaxTree();
            TaxNodeSelection selection = taxNodesGoal.get();
            selectedTaxNodes = selection.getSelected();
            excludedTaxNodes = selection.getExcluded();
            offPath.set(0);
            for (SmallTaxIdNode node : taxTree) {
                if (Rank.GENOME.equals(node.getRank())) {
                    map.put(node.getTaxId(), new long[pathLength(node)]);
                }
            }
            countUnstored = true;
            if (map.isEmpty()) {
                // Otherwise the goal reads every genome only to hand back an empty map without a word.
                // Without genome nodes no genome can be told from its taxon, so none counts as unstored
                // either: every one would, the stored ones included.
                if (getLogger().isWarnEnabled()) {
                    getLogger().warn("The database holds no GENOME-rank nodes - was it built with genomeNodes=true?");
                }
                countUnstored = false;
            }
            int stored = map.size();
            readFastas();
            if (getLogger().isInfoEnabled()) {
                getLogger().info("Path counts of " + stored + " stored genomes" + (countUnstored
                        ? " and " + (map.size() - stored) + " genomes the database did not store." : "."));
            }
            if (offPath.get() > 0 && getLogger().isWarnEnabled()) {
                getLogger().warn(offPath.get() + " k-mer occurrences of unstored genomes are held below the node the"
                        + " genome sits at and were left out.");
            }
            set(map);
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            map = null;
            store = null;
            taxTree = null;
            selectedTaxNodes = null;
            excludedTaxNodes = null;
        }
    }

    private static int pathLength(SmallTaxIdNode node) {
        int pathLength = 0;
        for (SmallTaxIdNode a = node; a != null; a = a.getParent()) {
            pathLength++;
        }
        return pathLength;
    }

    /**
     * Decides whether a contig is one the database's LCA update read, as {@code DBGoal.isContigInScope}
     * does: the scope restricts the RefSeq release only, and a fasta the project supplies itself always
     * takes part.
     *
     * @param node              the node the contig resolved to
     * @param fromRefSeqRelease whether the contig stems from the RefSeq release
     * @return whether the update read the contig
     */
    boolean isContigInScope(TaxIdNode node, boolean fromRefSeqRelease) {
        if (!fromRefSeqRelease) {
            return true;
        }
        switch (updateScope) {
            case OWN_TAXA_ONLY:
                return node != null && (selectedTaxNodes.isEmpty() || selectedTaxNodes.contains(node));
            case OTHER_TAXA_ONLY:
                return !selectedTaxNodes.contains(node);
            case ALL_BUT_EXCLUDED:
                return !excludedTaxNodes.contains(node);
            default:
                return true;
        }
    }

    @Override
    protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid) {
        // Counting the unstored genomes, every contig is a candidate and isContigInScope() decides;
        // otherwise only the selected taxa can hold a stored genome. An empty set means no restriction.
        return new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
                countUnstored ? Collections.emptySet() : selectedTaxNodes, selectedTaxNodes,
                isIncludeRefSeqFna() ? accessionMapGoal.get() : null, intConfigValue(GSConfigKey.KMER_SIZE),
                intConfigValue(GSConfigKey.MAX_DUST), intConfigValue(GSConfigKey.KMER_SAMPLING),
                booleanConfigValue(GSConfigKey.ASSEMBLY_ACCESSIONS_ONLY), contigsPerTaxid,
                booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES), (Rank) configValue(GSConfigKey.FOLD_TAXA_BELOW));
    }

    /**
     * FASTA reader that, for each k-mer of a contig, looks the k-mer up in the database and tallies it
     * against the current contig's genome, at the position along the genome's path the stored node
     * sits at.
     */
    protected class MyFastaReader extends AbstractStoreFastaReader {
        // Resolved once per contig by infoLine(), not once per k-mer - resolving it is a short walk of
        // the database's compact tree, and a contig's genome does not change while it is being read.
        // The node the genome sits at, index 0 of its counts, or null if the contig is not counted.
        private SmallTaxIdNode pathStart;
        // The genome's counts; for an unstored genome null until its first k-mer is found.
        private long[] genomeCounts;
        // The key of an unstored genome, or null for a stored one.
        private String unstoredKey;
        // The taxa selected for the database, which decide where the build folded a genome to.
        private final Set<TaxIdNode> selected;
        // The rank the database build folded genomes up to, or null - see foldUp().
        // Needed to find the node the build filed a genome under.
        private final Rank foldTaxaBelow;

        /**
         * Creates the reader.
         *
         * @param bufferSize             the read buffer size
         * @param taxNodes               the tax nodes whose contigs are read, empty for all
         * @param selected               the tax nodes selected for the database
         * @param accessionMap           the accession-to-taxid map
         * @param k                      the k-mer length
         * @param maxDust                the maximum allowed low-complexity (dust) run length
         * @param kMerSampling           the k-mer sampling step size
         * @param assemblyAccessionsOnly whether only genomic accessions are considered, dropping `NG_`, `NT_` and `NW_`
         * @param contigsPerTaxid        the per-taxid contig trie
         * @param enableLowerCaseBases   whether lower-case bases are included
         * @param foldTaxaBelow          the rank the database build folded genomes up to, or {@code null}
         */
        public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, Set<TaxIdNode> selected,
                             AccessionMap accessionMap, int k, int maxDust, int kMerSampling,
                             boolean assemblyAccessionsOnly, StringLong2DigitTrie contigsPerTaxid,
                             boolean enableLowerCaseBases, Rank foldTaxaBelow) {
            super(bufferSize, taxNodes, accessionMap, k, maxDust, kMerSampling, assemblyAccessionsOnly,
                    contigsPerTaxid, enableLowerCaseBases);
            this.selected = selected;
            this.foldTaxaBelow = foldTaxaBelow;
        }

        /**
         * Resolves and caches this contig's genome, in addition to the standard per-contig resolution
         * the superclass already does.
         */
        @Override
        protected void infoLine() {
            super.infoLine();
            if (includeContig && countUnstored && !isContigInScope(node, isRefSeqReleaseContig())) {
                includeContig = false;
            }
            updateGenomeNode();
        }

        /**
         * Resolves the node the current contig's genome sits at, one of the three cases of the class
         * comment. Starting from the tax node {@link #infoLine()} already resolved and folded as the
         * build folded it: if the database's tree has that node, the genome is looked for beneath it,
         * mirroring the genome-key cut {@code ReworkingStoreFastaReader} used when the database was
         * built, and if it is not there it sits at the node or its {@code DATA} child. If the tree does
         * not have the node, the genome sits at the nearest ancestor it does have. Leaves
         * {@link #pathStart} {@code null} when the contig is not counted.
         */
        private void updateGenomeNode() {
            pathStart = null;
            genomeCounts = null;
            unstoredKey = null;
            if (!includeContig || node == null || node.getTaxId() == null) {
                return;
            }
            // The build filed the genome under the folded node, not the one the accession names: with
            // foldTaxaBelow set, the latter is typically not even in the database's tree.
            TaxIdNode filedAt = ReworkingStoreFastaReader.foldUp(node, foldTaxaBelow, selected);
            SmallTaxIdNode n = taxTree.getNodeByTaxId(filedAt.getTaxId());
            int keyEnd = 1 + GenomeKeyTrie.genomeKeyLength(target, 1, accessionEnd());
            if (n != null) {
                SmallTaxIdNode genome = n;
                if (Rank.GENOME.ordinal() != genome.getRankOrdinal()) {
                    genome = n.getDescendantWithName(target, 1, keyEnd);
                }
                if (genome != null && Rank.GENOME.ordinal() == genome.getRankOrdinal()) {
                    // Stored: the genome's own node.
                    pathStart = genome;
                    genomeCounts = map.get(genome.getTaxId());
                    return;
                }
                if (!countUnstored) {
                    return;
                }
                // Named but not stored. The update merged the contig at the taxon's DATA child where
                // there is one, so that is where the k-mers it shares with the stored genomes went.
                SmallTaxIdNode data = Rank.DATA.ordinal() == n.getRankOrdinal() ? null : n.getDataChild();
                pathStart = data != null ? data : n;
            } else {
                if (!countUnstored) {
                    return;
                }
                // Not in the tree at all. Every node a k-mer of it can be stored at is a common ancestor
                // with something the tree has, so it lies at or above the nearest ancestor the tree has.
                for (TaxIdNode a = filedAt.getParent(); a != null && pathStart == null; a = a.getParent()) {
                    pathStart = taxTree.getNodeByTaxId(a.getTaxId());
                }
                if (pathStart == null) {
                    return;
                }
            }
            unstoredKey = unstoredGenomeKey(pathStart.getTaxId(),
                    new String(target, 1, keyEnd - 1, StandardCharsets.US_ASCII));
        }

        /**
         * Returns the exclusive end index of the accession within {@link #target}, the same cut
         * {@code ReworkingStoreFastaReader} uses to derive a genome's key from its accession.
         *
         * @return the exclusive end index of the accession
         */
        private int accessionEnd() {
            int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
            if (pos < 0) {
                pos = size;
                while (pos > 0 && (target[pos - 1] == '\n' || target[pos - 1] == '\r')) {
                    pos--;
                }
            }
            return pos;
        }

        /**
         * Looks the k-mer up in the database and, if found, tallies it against the current contig's
         * genome at the position on its path the stored node sits at.
         *
         * @return whether the k-mer was tallied
         */
        @Override
        protected boolean handleStore(long kmer) {
            if (pathStart == null || (genomeCounts == null && unstoredKey == null)) {
                // No genome resolved for this contig - nothing to tally against.
                return false;
            }
            String storedTaxId = store.getLong(kmer, null);
            if (storedTaxId == null) {
                // K-mer is not in the database - it did not land anywhere.
                return false;
            }
            SmallTaxIdNode storedNode = taxTree.getNodeByTaxId(storedTaxId);
            if (storedNode == null) {
                return false;
            }
            int index = 0;
            SmallTaxIdNode n = pathStart;
            while (n != null && n != storedNode) {
                n = n.getParent();
                index++;
            }
            if (n == null) {
                if (unstoredKey != null) {
                    // Below the node the genome sits at: a refined node of the finer tree, say. No row
                    // expects it there.
                    offPath.incrementAndGet();
                } else if (getLogger().isWarnEnabled()) {
                    // storedNode is not on the stored genome's own path to the root - should not happen
                    // given the lowest-common-ancestor update the database was built with, but guarded
                    // rather than indexed out of bounds.
                    getLogger().warn("K-mer stored at " + storedTaxId + " is not on the path of genome "
                            + pathStart.getTaxId() + ".");
                }
                return false;
            }
            long[] counts = genomeCounts;
            if (counts == null) {
                // The first k-mer of an unstored genome: its entry is made now, atomically, since another
                // reader may meet another contig of the same genome at the same time.
                SmallTaxIdNode start = pathStart;
                counts = map.computeIfAbsent(unstoredKey, key -> new long[pathLength(start)]);
                genomeCounts = counts;
            }
            // The array is shared by every reader that meets this genome - contigs of one genome may be
            // read by several threads at once - so locking on it yields per-genome mutual exclusion.
            synchronized (counts) {
                counts[index]++;
            }
            return true;
        }
    }
}
