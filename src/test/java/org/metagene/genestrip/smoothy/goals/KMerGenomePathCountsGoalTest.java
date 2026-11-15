package org.metagene.genestrip.smoothy.goals;

import org.junit.Assume;
import org.junit.BeforeClass;
import org.junit.Test;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.smoothy.SmoothyConfigKey;
import org.metagene.genestrip.smoothy.SmoothyGoalKey;
import org.metagene.genestrip.smoothy.SmoothyMaker;
import org.metagene.genestrip.smoothy.SmoothyProject;
import org.metagene.genestrip.smoothy.model.Topic;
import org.metagene.genestrip.smoothy.model.TopicModelResult;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.*;

/**
 * End-to-end tests of {@link KMerGenomePathCountsGoal}: a database is built from two synthetic genomes
 * with {@code genomeNodes=true}, and the goal's counts are compared with what a brute-force walk over
 * the same sequences predicts.
 * <p>
 * The genomes are random sequence, so none of their k-mers occurs in the RefSeq release the build
 * updates against, and where each k-mer ends up is decided by the two genomes alone: a k-mer of only
 * one genome stays at that genome's node (index {@code 0}), a k-mer of both goes to their lowest
 * common ancestor. Genome 1 is a two-contig WGS assembly carrying a repeat, so the tests also cover
 * that the contigs of one assembly add up under one genome and that a repeated k-mer is counted once
 * per occurrence.
 * <p>
 * The taxonomy and the RefSeq release are taken from the repository's {@code data/common}, since a
 * database build needs them and downloading them for a test would cost gigabytes; the tests are
 * skipped if they are not there.
 */
public class KMerGenomePathCountsGoalTest {
    private static final int K = 31;

    /**
     * Two contigs of one WGS assembly: their genome key is {@code NZ_ABCD}.
     */
    private static final String G1_CONTIG1 = "NZ_ABCD01000001.1";
    private static final String G1_CONTIG2 = "NZ_ABCD01000002.1";
    private static final String G1_KEY = "NZ_ABCD";
    /**
     * A finished replicon, which names its own genome, minus the version.
     */
    private static final String G2_CONTIG = "NC_999002.1";
    private static final String G2_KEY = "NC_999002";

    /**
     * Dengue virus 1 and 2, both "no rank" below 12637, which is below the species 3052464.
     */
    private static final String DENV1 = "11053";
    private static final String DENV2 = "11060";
    /**
     * Strains of the two, used to exercise {@code foldTaxaBelow}.
     */
    private static final String DENV1_STRAIN = "11054";
    private static final String DENV2_STRAIN = "11061";
    private static final String DENGUE_SPECIES = "3052464";
    /**
     * Dengue virus, the parent of the four serotypes, and serotype 4 with its genome in the release.
     */
    private static final String DENGUE_VIRUS = "12637";
    private static final String DENV4 = "11070";
    private static final String DENV4_ACCESSION = "NC_002640.1";
    private static final String DENV4_KEY = "NC_002640";
    /**
     * Genome 3: a stretch of DENV4, stored at DENV1.
     */
    private static final String G3_CONTIG = "NC_999003.1";
    private static final String G3_KEY = "NC_999003";
    private static final int DENV4_STRETCH_START = 1000;
    private static final int DENV4_STRETCH_END = 1400;

    private static final String[] G1_CONTIGS;
    private static final String G2_SEQ;

    static {
        Random random = new Random(4711);
        String shared = randomBases(random, 200);
        String repeat = randomBases(random, 80);
        G1_CONTIGS = new String[]{randomBases(random, 300) + repeat + shared, repeat + randomBases(random, 150)};
        G2_SEQ = randomBases(random, 250) + shared;
    }

    @BeforeClass
    public static void clearProjects() throws IOException {
        Assume.assumeTrue("Needs the taxonomy and the viral RefSeq release in data/common",
                new File(getCommonDir(), "nodes.dmp").exists()
                        && new File(getCommonDir(), "refseq/viral.1.1.genomic.fna.gz").exists());
        // Fresh databases on every run: one left over from an earlier version of the goal or of the
        // fixture would otherwise be loaded instead of built by the code under test.
        deleteRecursively(new File(getBaseDir(), "projects/smoothytest"));
        deleteRecursively(new File(getBaseDir(), "projects/smoothyfold"));
        deleteRecursively(new File(getBaseDir(), "projects/smoothyother"));
    }

    /**
     * Every k-mer occurrence is tallied at the right index: unique ones at the genome, shared ones at
     * the lowest common ancestor of the two genomes, repeats once per occurrence. The CSV file of
     * {@link KMerGenomePathCountsCSVGoal} is checked against the same counts here, since it needs the
     * same database and a maker of its own would cost another pass over the RefSeq catalog.
     */
    @Test
    public void countsMatchBruteForce() throws IOException {
        SmoothyMaker<SmoothyProject> maker = new SmoothyMaker<>(createProject("smoothytest", DENV1, DENV2, null, 0));
        try {
            SmallTaxTree tree = getDatabase(maker).getTaxTree();
            SmallTaxIdNode g1 = genomeNode(tree, DENV1, G1_KEY);
            SmallTaxIdNode g2 = genomeNode(tree, DENV2, G2_KEY);
            Map<String, long[]> counts = getCounts(maker);

            assertExpectedCounts(tree, g1, g2, counts);
            assertCsvMatches(maker, tree, counts);
            assertSerMatches(maker, counts);
        } finally {
            maker.dumpAll();
        }
    }

    /**
     * The CSV file holds one row per genome and path position, with the goal's count, the node at that
     * position, and a share that is the count over the genome's total.
     */
    private static void assertCsvMatches(SmoothyMaker<SmoothyProject> maker, SmallTaxTree tree,
                                         Map<String, long[]> counts) throws IOException {
        KMerGenomePathCountsCSVGoal<SmoothyProject> csvGoal = (KMerGenomePathCountsCSVGoal<SmoothyProject>) maker
                .getGoal(SmoothyGoalKey.KMER_GENOME_PATH_COUNTS_CSV);
        csvGoal.cleanThis();
        csvGoal.make();
        List<String> lines = Files.readAllLines(csvGoal.getFiles().get(0).toPath(), StandardCharsets.UTF_8);

        assertEquals("genome taxid;genome name;position;node taxid;node name;node rank;kmers;share;stored;",
                lines.get(0));
        int expectedRows = 0;
        for (long[] c : counts.values()) {
            expectedRows += c.length;
        }
        assertEquals(expectedRows, lines.size() - 1);

        Map<String, Double> shareSums = new HashMap<>();
        for (String line : lines.subList(1, lines.size())) {
            String[] cols = line.split(";", -1);
            String genomeTaxId = cols[0];
            int position = Integer.parseInt(cols[2]);
            long[] genomeCounts = counts.get(genomeTaxId);
            assertNotNull("Unknown genome " + genomeTaxId, genomeCounts);
            boolean unstored = KMerGenomePathCountsGoal.isUnstored(genomeTaxId);
            SmallTaxIdNode start = tree.getNodeByTaxId(KMerGenomePathCountsGoal.nodeOf(genomeTaxId));
            assertEquals(unstored ? KMerGenomePathCountsGoal.genomeKeyOf(genomeTaxId) : start.getName(), cols[1]);
            assertEquals(pathToRoot(start).get(position).getTaxId(), cols[3]);
            assertEquals(genomeCounts[position], Long.parseLong(cols[6]));
            assertEquals(String.valueOf(!unstored), cols[8]);
            if (!cols[7].isEmpty()) {
                shareSums.merge(genomeTaxId, Double.parseDouble(cols[7]), Double::sum);
            }
        }
        // Every genome with k-mers in the database has a distribution summing to one over its path.
        for (Map.Entry<String, Double> e : shareSums.entrySet()) {
            assertEquals(e.getKey(), 1.0, e.getValue(), 1e-6);
        }
        assertTrue(shareSums.size() >= 2);
    }

    /**
     * The SER file reads back as exactly the goal's counts.
     */
    private static void assertSerMatches(SmoothyMaker<SmoothyProject> maker, Map<String, long[]> counts)
            throws IOException {
        KMerGenomePathCountsSERGoal<SmoothyProject> serGoal = (KMerGenomePathCountsSERGoal<SmoothyProject>) maker
                .getGoal(SmoothyGoalKey.KMER_GENOME_PATH_COUNTS_SER);
        serGoal.cleanThis();
        serGoal.make();
        KMerGenomePathCountsSERGoal.StoredPathCounts stored;
        try {
            stored = KMerGenomePathCountsSERGoal.load(serGoal.getFiles().get(0));
        } catch (ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
        // The file names the database it was measured on.
        assertNotNull(stored.getDbMD5());
        assertEquals(getDatabase(maker).getConfigInfo().getProperty(GSProject.DB_MD5), stored.getDbMD5());
        Map<String, long[]> loaded = stored.getCounts();
        assertEquals(counts.keySet(), loaded.keySet());
        for (String taxId : counts.keySet()) {
            assertArrayEquals(taxId, counts.get(taxId), loaded.get(taxId));
        }
    }

    /**
     * The reader threads share the counters, so reading in parallel must not change a single count.
     */
    @Test
    public void multiThreadedEqualsSingleThreaded() throws IOException {
        Map<String, long[]> single;
        SmoothyMaker<SmoothyProject> maker = new SmoothyMaker<>(createProject("smoothytest", DENV1, DENV2, null, 0));
        try {
            single = getCounts(maker);
        } finally {
            maker.dumpAll();
        }
        maker = new SmoothyMaker<>(createProject("smoothytest", DENV1, DENV2, null, 4));
        try {
            Map<String, long[]> multi = getCounts(maker);
            assertEquals(single.keySet(), multi.keySet());
            for (String taxId : single.keySet()) {
                assertArrayEquals(taxId, single.get(taxId), multi.get(taxId));
            }
        } finally {
            maker.dumpAll();
        }
    }

    /**
     * There is one entry per genome node of the database, and each array spans the genome's whole
     * path to the root - so that every node a k-mer can be stored at has a slot. Entries of genomes the
     * database did not store come on top, each spanning the path from the node it sits at.
     */
    @Test
    public void oneEntryPerGenomeSpanningItsPath() throws IOException {
        SmoothyMaker<SmoothyProject> maker = new SmoothyMaker<>(createProject("smoothytest", DENV1, DENV2, null, 0));
        try {
            SmallTaxTree tree = getDatabase(maker).getTaxTree();
            Map<String, long[]> result = getCounts(maker);
            Set<String> genomes = new HashSet<>();
            for (SmallTaxIdNode node : tree) {
                if (Rank.GENOME.equals(node.getRank())) {
                    genomes.add(node.getTaxId());
                    assertEquals(pathToRoot(node).size(), result.get(node.getTaxId()).length);
                }
            }
            for (String key : result.keySet()) {
                if (KMerGenomePathCountsGoal.isUnstored(key)) {
                    SmallTaxIdNode start = tree.getNodeByTaxId(KMerGenomePathCountsGoal.nodeOf(key));
                    assertNotNull(key, start);
                    assertNotSame(key, Rank.GENOME, start.getRank());
                    assertEquals(key, pathToRoot(start).size(), result.get(key).length);
                    // Only genomes with something in the database get an entry.
                    assertTrue(key, Arrays.stream(result.get(key)).sum() > 0);
                }
            }
            assertEquals(genomes, storedKeys(result));
            // Both contigs of the WGS assembly went into one genome, and nothing was filed per contig.
            SmallTaxIdNode denv1 = tree.getNodeByTaxId(DENV1);
            assertTrue(genomes.contains(genomeNode(tree, DENV1, G1_KEY).getTaxId()));
            assertNull(denv1.getDescendantWithName(G1_CONTIG1.substring(0, G1_CONTIG1.indexOf('.'))));
            assertNull(denv1.getDescendantWithName(G1_CONTIG2.substring(0, G1_CONTIG2.indexOf('.'))));
        } finally {
            maker.dumpAll();
        }
    }

    /**
     * With {@code foldTaxaBelow} the build files a strain's genome under its species - a node the
     * accession map does not name, while the strain itself is not even in the database's tree. The
     * goal must look for the genome where the build put it.
     */
    @Test
    public void genomesFoldedToTheSpeciesAreFound() throws IOException {
        SmoothyMaker<SmoothyProject> maker = new SmoothyMaker<>(
                createProject("smoothyfold", DENV1_STRAIN, DENV2_STRAIN, Rank.SPECIES, 0));
        try {
            SmallTaxTree tree = getDatabase(maker).getTaxTree();
            SmallTaxIdNode g1 = genomeNode(tree, DENGUE_SPECIES, G1_KEY);
            SmallTaxIdNode g2 = genomeNode(tree, DENGUE_SPECIES, G2_KEY);
            assertEquals(DENGUE_SPECIES, g1.getParent().getTaxId());
            assertEquals(DENGUE_SPECIES, g2.getParent().getTaxId());
            Map<String, long[]> result = getCounts(maker);

            assertExpectedCounts(tree, g1, g2, result);
            // Shared k-mers sit at the species, one step above each genome.
            assertTrue(result.get(g1.getTaxId())[1] > 0);
        } finally {
            maker.dumpAll();
        }
    }

    /**
     * A new maker - as on the next run from the command line - gets the counts from the SER file,
     * without counting again, and gets exactly what was counted.
     */
    @Test
    public void loadGoalReadsTheSerFileWithoutCounting() throws IOException {
        Map<String, long[]> counted;
        SmoothyMaker<SmoothyProject> maker = new SmoothyMaker<>(createProject("smoothytest", DENV1, DENV2, null, 0));
        try {
            FileGoal<SmoothyProject> serGoal = (FileGoal<SmoothyProject>) maker.getGoal(SmoothyGoalKey.KMER_GENOME_PATH_COUNTS_SER);
            serGoal.cleanThis();
            // Taken before the SER file is written: once it is, its only dependent is made and
            // Genestrip frees the counts, so asking for them afterwards would count all over again.
            counted = getCounts(maker);
            serGoal.make();
        } finally {
            maker.dumpAll();
        }

        maker = new SmoothyMaker<>(createProject("smoothytest", DENV1, DENV2, null, 0));
        try {
            Map<String, long[]> loaded = getLoadedCounts(maker);
            assertFalse("The load goal must not count again",
                    maker.getGoal(SmoothyGoalKey.KMER_GENOME_PATH_COUNTS).isMade());
            assertEquals(counted.keySet(), loaded.keySet());
            for (String taxId : counted.keySet()) {
                assertArrayEquals(taxId, counted.get(taxId), loaded.get(taxId));
            }
        } finally {
            maker.dumpAll();
        }
    }

    /**
     * Counts recorded against another database - one rebuilt since - are refused, not used.
     */
    @Test
    public void loadGoalRefusesCountsOfAnotherDatabase() throws IOException {
        SmoothyMaker<SmoothyProject> maker = new SmoothyMaker<>(createProject("smoothytest", DENV1, DENV2, null, 0));
        File serFile = ((FileGoal<SmoothyProject>) maker.getGoal(SmoothyGoalKey.KMER_GENOME_PATH_COUNTS_SER)).getFile();
        String otherMD5 = "00000000000000000000000000000000";
        try {
            // The database must exist for the check to have something to compare against.
            String dbMD5 = getDatabase(maker).getConfigInfo().getProperty(GSProject.DB_MD5);
            serFile.getParentFile().mkdirs();
            try (ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(serFile))) {
                out.writeObject(otherMD5);
                out.writeObject(new HashMap<String, long[]>());
            }
            try {
                getLoadedCounts(maker);
                fail("Counts of another database must be refused");
            } catch (StalePathCountsException e) {
                // Checked by what was compared, not by the wording of the message.
                assertEquals(serFile, e.getCountsFile());
                assertEquals(otherMD5, e.getCountsMD5());
                assertEquals(dbMD5, e.getDbMD5());
            }
        } finally {
            maker.dumpAll();
            // Not to leave a file behind that every later run would refuse.
            serFile.delete();
        }
    }

    /**
     * The whole chain on a sample of reads from genome 1 alone: match, path counts, topic model, CSV.
     * Genome 1 must come out as the topic carrying the sample, and genome 2, of which the sample holds
     * nothing of its own, must not be a topic at all.
     */
    @Test
    public void topicModelFindsTheGenomeASampleComesFrom() throws IOException {
        // Reads of 100 bases every 5 positions of both contigs of genome 1, and nothing else.
        StringBuilder fastq = new StringBuilder();
        int read = 0;
        for (String contig : G1_CONTIGS) {
            for (int i = 0; i + 100 <= contig.length(); i += 5) {
                fastq.append("@read").append(read++).append('\n').append(contig, i, i + 100).append("\n+\n");
                for (int j = 0; j < 100; j++) {
                    fastq.append('I');
                }
                fastq.append('\n');
            }
        }
        File sample = new File(getBaseDir(), "smoothy_g1_sample.fastq");
        write(sample, fastq.toString());

        SmoothyMaker<SmoothyProject> maker = new SmoothyMaker<>(createProject("smoothytest", DENV1, DENV2, null, 0,
                "g1sample", new String[]{sample.getAbsolutePath()}));
        try {
            SmallTaxTree tree = getDatabase(maker).getTaxTree();
            String g1 = genomeNode(tree, DENV1, G1_KEY).getTaxId();
            String g2 = genomeNode(tree, DENV2, G2_KEY).getTaxId();

            Map<String, TopicModelResult> results = ((ObjectGoal<Map<String, TopicModelResult>, SmoothyProject>) maker
                    .getGoal(SmoothyGoalKey.TOPIC_MODEL)).get();
            assertEquals(1, results.size());
            TopicModelResult result = results.values().iterator().next();

            Topic top = result.getTopics().get(0);
            assertEquals(g1, top.getNode());
            assertEquals(1, top.getGenomes());
            assertTrue(top.getKmerShare() > 0.9);
            for (Topic topic : result.getTopics()) {
                assertNotEquals(g2, topic.getNode());
            }
            assertTrue(result.getUncoveredNodes().isEmpty());

            FileListGoal<SmoothyProject> csvGoal = (FileListGoal<SmoothyProject>) maker.getGoal(SmoothyGoalKey.TOPIC_MODEL_CSV);
            csvGoal.cleanThis();
            csvGoal.make();
            List<String> lines = Files.readAllLines(csvGoal.getFiles().get(0).toPath(), StandardCharsets.UTF_8);
            assertEquals("taxid;name;rank;kmers;genomes;kmer share;expected hold;abundance;consistency weight;", lines.get(0));
            assertTrue(lines.get(1), lines.get(1).startsWith(g1 + ";" + G1_KEY + ";GENOME;"));
            assertEquals(result.getTopics().size() + 1, lines.size());
        } finally {
            maker.dumpAll();
        }
    }

    /**
     * An organism outside the database, the second case of the paper's p. 21: the database holds Dengue
     * virus 1 and 2 only, and Dengue virus 4 of the RefSeq release is not in its tree. Genome 3, stored
     * at DENV1, carries a stretch of DENV4, so the database holds k-mers of DENV4 - at a common ancestor,
     * since DENV4 took part in the update.
     * <p>
     * DENV4 must then get an entry of its own, sitting at the nearest node the tree has above its taxon,
     * and hold at least what genome 3 holds of the shared stretch at the node it is stored at. Counting
     * the stored genomes only must leave their counts exactly as they are and write files of another
     * name.
     */
    @Test
    public void genomeOutsideTheDatabaseIsCountedAtTheNearestAncestor() throws IOException {
        Map<String, long[]> withUnstored;
        SmoothyMaker<SmoothyProject> maker = new SmoothyMaker<>(createOtherProject(true, null, null));
        try {
            SmallTaxTree tree = getDatabase(maker).getTaxTree();
            assertNull("DENV4 must not be in the database", tree.getNodeByTaxId(DENV4));
            SmallTaxIdNode g3 = genomeNode(tree, DENV1, G3_KEY);
            withUnstored = getCounts(maker);

            String denv4 = unstoredKey(withUnstored, DENV4_KEY);
            SmallTaxIdNode start = tree.getNodeByTaxId(KMerGenomePathCountsGoal.nodeOf(denv4));
            // Dengue virus, the parent of all four serotypes and the nearest node the tree has.
            assertEquals(DENGUE_VIRUS, start.getTaxId());
            assertEquals(pathToRoot(start).size(), withUnstored.get(denv4).length);

            // Where genome 3 holds the stretch it shares with DENV4: above DENV1, on DENV4's path too.
            long[] g3Counts = withUnstored.get(g3.getTaxId());
            int g3Index = 0;
            for (int i = 1; i < g3Counts.length; i++) {
                if (g3Counts[i] > g3Counts[g3Index]) {
                    g3Index = i;
                }
            }
            SmallTaxIdNode shared = pathToRoot(g3).get(g3Index);
            assertTrue("The stretch of DENV4 must not stay at genome 3", g3Index > 0);
            assertTrue(g3Counts[g3Index] > (DENV4_STRETCH_END - DENV4_STRETCH_START) / 2);
            int denv4Index = pathToRoot(start).indexOf(shared);
            assertTrue(denv4Index >= 0);
            assertTrue(withUnstored.get(denv4)[denv4Index] >= g3Counts[g3Index]);

            FileGoal<SmoothyProject> csvGoal = (FileGoal<SmoothyProject>) maker.getGoal(SmoothyGoalKey.KMER_GENOME_PATH_COUNTS_CSV);
            assertTrue(csvGoal.getFiles().get(0).getName().contains("_withunstored"));
            assertCsvMatches(maker, tree, withUnstored);
        } finally {
            maker.dumpAll();
        }

        maker = new SmoothyMaker<>(createOtherProject(false, null, null));
        try {
            Map<String, long[]> storedOnly = getCounts(maker);
            assertEquals(storedKeys(withUnstored), storedOnly.keySet());
            for (String taxId : storedOnly.keySet()) {
                assertArrayEquals(taxId, withUnstored.get(taxId), storedOnly.get(taxId));
            }
            FileGoal<SmoothyProject> serGoal = (FileGoal<SmoothyProject>) maker.getGoal(SmoothyGoalKey.KMER_GENOME_PATH_COUNTS_SER);
            assertFalse(serGoal.getFiles().get(0).getName().contains("_withunstored"));
        } finally {
            maker.dumpAll();
        }
    }

    /**
     * The mixed case: a sample of DENV4, which the database does not hold, beside genome 1, which it
     * does. Counting DENV4 puts it into the pool of the topic at Dengue virus, which therefore speaks for
     * one genome more than with the stored genomes alone, and that topic carries a real share of the
     * sample.
     */
    @Test
    public void genomeOutsideTheDatabaseJoinsThePoolOfItsAncestorsTopic() throws IOException {
        StringBuilder fastq = new StringBuilder();
        int read = appendReads(fastq, 0, releaseSequence(DENV4_ACCESSION));
        for (String contig : G1_CONTIGS) {
            read = appendReads(fastq, read, contig);
        }
        File sample = new File(getBaseDir(), "smoothy_denv4_g1_sample.fastq");
        write(sample, fastq.toString());
        String[] fastqFiles = {sample.getAbsolutePath()};

        Topic withUnstored;
        SmoothyMaker<SmoothyProject> maker = new SmoothyMaker<>(createOtherProject(true, "denv4g1", fastqFiles));
        try {
            withUnstored = topic(topicModel(maker), DENGUE_VIRUS);
        } finally {
            maker.dumpAll();
        }
        Topic storedOnly;
        maker = new SmoothyMaker<>(createOtherProject(false, "denv4g1", fastqFiles));
        try {
            storedOnly = topic(topicModel(maker), DENGUE_VIRUS);
        } finally {
            maker.dumpAll();
        }
        assertNotNull("Dengue virus must be a topic", withUnstored);
        assertTrue(withUnstored.getGenomes() > (storedOnly == null ? 0 : storedOnly.getGenomes()));
        assertTrue(withUnstored.getKmerShare() > 0.1);
    }

    private static TopicModelResult topicModel(SmoothyMaker<SmoothyProject> maker) {
        Map<String, TopicModelResult> results = ((ObjectGoal<Map<String, TopicModelResult>, SmoothyProject>) maker
                .getGoal(SmoothyGoalKey.TOPIC_MODEL)).get();
        assertEquals(1, results.size());
        return results.values().iterator().next();
    }

    private static Topic topic(TopicModelResult result, String node) {
        for (Topic t : result.getTopics()) {
            if (t.getNode().equals(node)) {
                return t;
            }
        }
        return null;
    }

    /**
     * Appends reads of 100 bases every 5 positions of the sequence.
     *
     * @return the number of the next read
     */
    private static int appendReads(StringBuilder fastq, int read, String seq) {
        for (int i = 0; i + 100 <= seq.length(); i += 5) {
            fastq.append("@read").append(read++).append('\n').append(seq, i, i + 100).append("\n+\n");
            for (int j = 0; j < 100; j++) {
                fastq.append('I');
            }
            fastq.append('\n');
        }
        return read;
    }

    private static Set<String> storedKeys(Map<String, long[]> counts) {
        Set<String> res = new HashSet<>();
        for (String key : counts.keySet()) {
            if (!KMerGenomePathCountsGoal.isUnstored(key)) {
                res.add(key);
            }
        }
        return res;
    }

    private static String unstoredKey(Map<String, long[]> counts, String genomeKey) {
        String res = null;
        for (String key : counts.keySet()) {
            if (genomeKey.equals(KMerGenomePathCountsGoal.genomeKeyOf(key))) {
                assertNull("Two entries for " + genomeKey, res);
                res = key;
            }
        }
        assertNotNull("No entry for the unstored genome " + genomeKey, res);
        return res;
    }

    /**
     * The sequence of a contig of the viral RefSeq release, upper-cased.
     */
    private static String releaseSequence(String accession) throws IOException {
        StringBuilder seq = new StringBuilder();
        boolean in = false;
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(
                new FileInputStream(new File(getCommonDir(), "refseq/viral.1.1.genomic.fna.gz"))),
                StandardCharsets.US_ASCII))) {
            for (String line = reader.readLine(); line != null; line = reader.readLine()) {
                if (line.startsWith(">")) {
                    if (in) {
                        break;
                    }
                    in = line.startsWith(">" + accession + " ");
                } else if (in) {
                    seq.append(line.trim().toUpperCase(Locale.ROOT));
                }
            }
        }
        assertTrue("No " + accession + " in the viral release", seq.length() > 0);
        return seq.toString();
    }

    /**
     * Sets up a project whose database is selected as Dengue virus 1 and 2 only, so that DENV3 and
     * DENV4 of the release are outside its tree. It holds genomes 1 and 2 as the other projects do, and
     * genome 3 at DENV1, which is a stretch of DENV4.
     *
     * @param unstoredGenomes {@code smoothy.unstoredGenomes}
     * @param sampleKey       the key the samples are analysed under, or {@code null}
     * @param fastqFiles      the samples, or {@code null}
     */
    private static SmoothyProject createOtherProject(boolean unstoredGenomes, String sampleKey, String[] fastqFiles)
            throws IOException {
        SmoothyProject project = createProject("smoothyother", DENV1, DENV2, null, 0, sampleKey, fastqFiles);
        write(new File(project.getProjectDir(), "taxids.txt"), DENV1 + "\n" + DENV2 + "\n");
        write(new File(project.getProjectDir(), "additional.txt"),
                DENV1 + " g1.fasta\n" + DENV2 + " g2.fasta\n" + DENV1 + " g3.fasta\n");
        write(new File(project.getFastaDir(), "g3.fasta"), ">" + G3_CONTIG + " synthetic\n"
                + releaseSequence(DENV4_ACCESSION).substring(DENV4_STRETCH_START, DENV4_STRETCH_END) + "\n");
        project.initConfigParam(SmoothyConfigKey.UNSTORED_GENOMES, unstoredGenomes);
        return project;
    }

    private static Map<String, long[]> getLoadedCounts(SmoothyMaker<SmoothyProject> maker) {
        return ((ObjectGoal<Map<String, long[]>, SmoothyProject>) maker
                .getGoal(SmoothyGoalKey.LOAD_KMER_GENOME_PATH_COUNTS)).get();
    }

    private static void assertExpectedCounts(SmallTaxTree tree, SmallTaxIdNode g1, SmallTaxIdNode g2,
                                             Map<String, long[]> result) {
        SmallTaxIdNode lca = tree.getLowestCommonAncestor(g1, g2);
        int lcaIndex1 = pathToRoot(g1).indexOf(lca);
        int lcaIndex2 = pathToRoot(g2).indexOf(lca);

        long[] expected1 = new long[pathToRoot(g1).size()];
        Set<String> g2Kmers = kmersOf(G2_SEQ);
        for (String contig : G1_CONTIGS) {
            tally(contig, g2Kmers, lcaIndex1, expected1);
        }
        long[] expected2 = new long[pathToRoot(g2).size()];
        tally(G2_SEQ, kmersOf(G1_CONTIGS), lcaIndex2, expected2);

        // Sanity check of the fixture itself: both slots are actually exercised.
        assertTrue(expected1[0] > 0 && expected1[lcaIndex1] > 0);
        assertArrayEquals(expected1, result.get(g1.getTaxId()));
        assertArrayEquals(expected2, result.get(g2.getTaxId()));
    }

    /**
     * Adds every k-mer occurrence of the sequence to {@code counts}: at {@code sharedIndex} if the other
     * genome has the k-mer too, at {@code 0} otherwise.
     */
    private static void tally(String seq, Set<String> otherKmers, int sharedIndex, long[] counts) {
        for (int i = 0; i + K <= seq.length(); i++) {
            counts[otherKmers.contains(canonical(seq.substring(i, i + K))) ? sharedIndex : 0]++;
        }
    }

    private static Set<String> kmersOf(String... seqs) {
        Set<String> res = new HashSet<>();
        for (String seq : seqs) {
            for (int i = 0; i + K <= seq.length(); i++) {
                res.add(canonical(seq.substring(i, i + K)));
            }
        }
        return res;
    }

    /**
     * The lexicographically smaller of the k-mer and its reverse complement.
     */
    private static String canonical(String kmer) {
        StringBuilder rc = new StringBuilder(kmer.length());
        for (int i = kmer.length() - 1; i >= 0; i--) {
            switch (kmer.charAt(i)) {
                case 'A':
                    rc.append('T');
                    break;
                case 'C':
                    rc.append('G');
                    break;
                case 'G':
                    rc.append('C');
                    break;
                default:
                    rc.append('A');
                    break;
            }
        }
        String r = rc.toString();
        return kmer.compareTo(r) <= 0 ? kmer : r;
    }

    private static String randomBases(Random random, int length) {
        char[] bases = {'A', 'C', 'G', 'T'};
        StringBuilder sb = new StringBuilder(length);
        for (int i = 0; i < length; i++) {
            sb.append(bases[random.nextInt(4)]);
        }
        return sb.toString();
    }

    /**
     * The node itself at index 0, then each ancestor up to the root - the goal's indexing.
     */
    private static List<SmallTaxIdNode> pathToRoot(SmallTaxIdNode node) {
        List<SmallTaxIdNode> path = new ArrayList<>();
        for (SmallTaxIdNode n = node; n != null; n = n.getParent()) {
            path.add(n);
        }
        return path;
    }

    private static SmallTaxIdNode genomeNode(SmallTaxTree tree, String underTaxId, String key) {
        SmallTaxIdNode parent = tree.getNodeByTaxId(underTaxId);
        assertNotNull("Tax id " + underTaxId + " is not in the database", parent);
        SmallTaxIdNode genome = parent.getDescendantWithName(key);
        assertNotNull("No genome node " + key + " below " + underTaxId, genome);
        assertSame(Rank.GENOME, genome.getRank());
        return genome;
    }

    private static Database getDatabase(SmoothyMaker<SmoothyProject> maker) {
        return ((ObjectGoal<Database, SmoothyProject>) maker.getGoal(GSGoalKey.LOAD_DB)).get();
    }

    private static Map<String, long[]> getCounts(SmoothyMaker<SmoothyProject> maker) {
        return ((ObjectGoal<Map<String, long[]>, SmoothyProject>) maker
                .getGoal(SmoothyGoalKey.KMER_GENOME_PATH_COUNTS)).get();
    }

    /**
     * Sets up a project whose database holds genome 1 filed at {@code taxId1} and genome 2 at
     * {@code taxId2}, both read as additional fastas, besides the Dengue genomes of the RefSeq release.
     * <p>
     * The release genomes cannot be kept out without {@code refseq.filldb=false}, which fetches the
     * same taxa from Genbank instead. They do no harm: none of their k-mers occurs in the random
     * genomes, so they change nothing the tests predict.
     */
    private static SmoothyProject createProject(String name, String taxId1, String taxId2, Rank foldTaxaBelow,
                                                int threads) throws IOException {
        return createProject(name, taxId1, taxId2, foldTaxaBelow, threads, null, null);
    }

    /**
     * Same as {@link #createProject(String, String, String, Rank, int)}, with samples to analyse.
     *
     * @param sampleKey  the key the samples are analysed under, as {@code -k} gives it
     * @param fastqFiles the samples, as {@code -f} gives them
     */
    private static SmoothyProject createProject(String name, String taxId1, String taxId2, Rank foldTaxaBelow,
                                                int threads, String sampleKey, String[] fastqFiles) throws IOException {
        GSCommon common = new GSCommon(getBaseDir()) {
            @Override
            public File getCommonDir() {
                return KMerGenomePathCountsGoalTest.getCommonDir();
            }
        };
        SmoothyProject project = new SmoothyProject(common, name, sampleKey, fastqFiles, null, null, null, null, null,
                null, null, true);
        File fastaDir = project.getFastaDir();
        write(new File(project.getProjectDir(), "taxids.txt"), DENGUE_SPECIES + "\n");
        write(new File(project.getProjectDir(), "categories.txt"), "viral\n");
        write(new File(project.getProjectDir(), "additional.txt"), taxId1 + " g1.fasta\n" + taxId2 + " g2.fasta\n");
        write(new File(fastaDir, "g1.fasta"), ">" + G1_CONTIG1 + " synthetic\n" + G1_CONTIGS[0] + "\n>"
                + G1_CONTIG2 + " synthetic\n" + G1_CONTIGS[1] + "\n");
        write(new File(fastaDir, "g2.fasta"), ">" + G2_CONTIG + " synthetic\n" + G2_SEQ + "\n");

        project.initConfigParam(GSConfigKey.GENOME_NODES, true);
        project.initConfigParam(GSConfigKey.THREADS, threads);
        project.initConfigParam(GSConfigKey.PROGRESS_BAR, false);
        if (foldTaxaBelow != null) {
            project.initConfigParam(GSConfigKey.FOLD_TAXA_BELOW, foldTaxaBelow);
        }
        return project;
    }

    private static void write(File file, String content) throws IOException {
        file.getParentFile().mkdirs();
        Files.write(file.toPath(), content.getBytes(StandardCharsets.UTF_8));
    }

    private static void deleteRecursively(File dir) throws IOException {
        if (dir.exists()) {
            try (Stream<Path> paths = Files.walk(dir.toPath())) {
                paths.sorted(Comparator.reverseOrder()).map(Path::toFile).forEach(File::delete);
            }
        }
    }

    /**
     * Below the build output, so that the test projects never land beside the sources.
     */
    private static File getBaseDir() {
        String buildDir = System.getProperty("project.build.directory");
        File target = buildDir != null ? new File(buildDir)
                : new File(KMerGenomePathCountsGoalTest.class.getProtectionDomain().getCodeSource().getLocation()
                .getFile()).getParentFile();
        return new File(target, "data");
    }

    /**
     * Genestrip's shared {@code data/common}, which holds the taxonomy and the RefSeq downloads: in
     * the data directory the {@code data.directory} system property names, otherwise in this
     * repository's {@code data}.
     */
    private static File getCommonDir() {
        String dataDir = System.getProperty("data.directory");
        File data = dataDir != null ? new File(dataDir)
                : new File(getBaseDir().getParentFile().getParentFile(), "data");
        return new File(data, "common");
    }
}
