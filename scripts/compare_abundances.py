#!/usr/bin/env python3
"""Compares abundance estimates of one simulated sample against its ground truth.

Two levels (`--level`):

* genome  (default) - genome nodes of a Genestrip database built with genomeNodes=true.
* species - every estimate and the truth are summed up to their species. This is the level Kraken2 and
            Bracken work on, so they can only be compared here. Needs `--taxonomy` and `--library`.

Methods compared:

* naive     - Genestrip's classified reads (`reads` column of the `match` CSV). On genome level only the reads
              at genome nodes count, reads Genestrip places higher in the tree (genus, family ...) cannot be
              given to a genome. On species level the reads at and below a species count for it, those
              above a species cannot be given to one. Two variants, both reported:
              - `naive renorm.`: the usable reads are rescaled to sum to one, so the dropped reads are spread
                proportionally. The kindest reading for the naive method.
              - `naive`: shares of all classified reads; the dropped reads count as error, exactly as the
                mass of taxon topics does for the topic model. This is the like-for-like comparison.
* topic models - one or more `topicmodelcsv` files (e.g. with and without consistency weight). Topics at
              taxa rather than genomes (on species level: above a species) carry mass nobody gets, which
              counts as error.
* Kraken2   - (species level, `--kreport`) the clade reads of each species in the Kraken2 report, the
              Kraken2 counterpart of `naive`; `Kraken2 renorm.` and `Kraken2` as above.
* Bracken   - (species level, `--bracken`) `new_est_reads` of the Bracken output. `Bracken renorm.` is
              Bracken's own `fraction_total_reads`, i.e. the species reads rescaled to sum to one, which is
              what Bracken reports to its users. `Bracken` takes the shares of all reads Kraken2 classified,
              so the reads Bracken could not push down to a species count as error, as for `naive`.

The Kraken2 and Bracken database must be built from the same sequences and taxonomy as the Genestrip
database (see kraken_library.py), else the comparison measures the databases, not the methods. Kraken2
counts a read pair once, Genestrip and the truth count reads; this does not matter, since only shares are
compared.

Mass a method holds but cannot give to any genome or species ("mass not at a genome/species") is added to its
L1, since the truth puts all mass there. For the naive organism share this also holds the reads of genomes
without known length (no fasta, no topic): they cannot be turned into organisms and count as error. The
dropped reads of the strict naive variant enter the organism share with their read share, as there is no
length to convert them by.

Genomes are identified the way Genestrip names its genome nodes: by accession for databases with one sequence
per genome (viral), by tax id for genomes of many sequences whose ids carry `|kraken:taxid|` (the
per-accession FASTA files of `extractrefseqcsv`, e.g. protozoa). Fastq and fasta files may be gzipped.

Two kinds of shares are compared on genome level, since they answer different questions:

* read share     - share of the sample's sequence a genome accounts for.
                   truth: reads per genome from the ISS read names; naive: read shares;
                   topic model: `kmer share` (theta).
* organism share - share of organisms (genome copies), i.e. read share divided by genome length.
                   truth: reads / length; naive: reads / length; topic model: `abundance` (a_m).

On species level only the read share is compared: Kraken2 and Bracken give no genome lengths, and a species
may hold several genomes of different length.

The ground truth is counted from the reads themselves: InSilicoSeq names every read after its source
sequence (`@NC_040360.1_0_0/1`), so no assumption about how ISS interprets its abundance or coverage file is
needed.

Distance per method: L1 = sum over all genomes (species) of |estimate - truth| (0 = perfect, 2 = disjoint).
Half of it is the share of mass assigned wrongly (total variation, equal to Bray-Curtis for compositions).

Examples:
    python compare_abundances.py --fastq ../data/projects/viral/fastq/miseq_cov_R1.fastq \
        --genomes ../data/projects/viral/insilicoseq/miseq_ncbi_ncbi_genomes.fasta \
        --match ../data/projects/viral/eval/weighted/viral_match_miseq_cov.csv \
        --topics weighted=../data/projects/viral/eval/weighted/viral_topicmodelcsv_miseq_cov.csv \
                 unweighted=../data/projects/viral/eval/unweighted/viral_topicmodelcsv_miseq_cov.csv \
        --out miseq_cov_comparison.csv

    python compare_abundances.py --level species --fastq ../data/projects/viral/fastq/miseq_cov_R1.fastq \
        --taxonomy ../data/common/nodes.dmp --library ../data/kraken/viral_library.fna \
        --match ../data/projects/viral/eval/weighted/viral_match_miseq_cov.csv \
        --topics weighted=../data/projects/viral/eval/weighted/viral_topicmodelcsv_miseq_cov.csv \
        --kreport ../data/projects/viral/kraken/miseq_cov.kreport \
        --bracken ../data/projects/viral/kraken/miseq_cov_S.bracken
"""
import argparse
import csv
import gzip
import re
import sys
from collections import Counter

READ_NAME = re.compile(r"^@(.+)_\d+_\d+(/[12])?$")
KRAKEN_TAXID = re.compile(r"\|kraken:taxid\|(\d+)")


def open_text(path, **kwargs):
    return gzip.open(path, "rt", **kwargs) if path.endswith(".gz") else open(path, **kwargs)


def accession(seq_id):
    """NC_040360.1 -> NC_040360, the name Genestrip gives a genome node."""
    return seq_id.split(".")[0]


def genome_of(seq_id):
    """The genome node a sequence belongs to.

    Two cases, told apart by the sequence id:
    * `NC_000521.4|kraken:taxid|36329` (from `extractrefseqcsv`, e.g. protozoa): a genome consists of many
      sequences, and Genestrip names its genome node after the tax id, so the genome is the tax id.
    * `NC_040360.1` (e.g. viral, one sequence per genome): Genestrip names the genome node after the
      accession without version.
    """
    m = KRAKEN_TAXID.search(seq_id)
    return m.group(1) if m else accession(seq_id)


def truth_reads(fastq):
    """Counts the reads per source genome in an ISS fastq file (one mate file is enough)."""
    counts = Counter()
    with open_text(fastq, encoding="ascii") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                m = READ_NAME.match(line.strip())
                if not m:
                    sys.exit("Not an InSilicoSeq read name: " + line.strip())
                counts[genome_of(m.group(1))] += 1
    return counts


def fasta_lengths(paths):
    """Length per genome, summed over all its sequences."""
    lengths = Counter()
    for path in paths:
        name = None
        with open_text(path, encoding="ascii") as f:
            for line in f:
                if line.startswith(">"):
                    name = genome_of(line[1:].split()[0])
                    lengths[name] += 0
                elif name is not None:
                    lengths[name] += len(line.strip())
    return dict(lengths)


def library_taxids(path):
    """Accession -> taxid from the `kraken:taxid` headers of a Kraken2 library fasta (kraken_library.py)."""
    taxids = {}
    with open_text(path, encoding="ascii", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                m = KRAKEN_TAXID.search(line)
                if m:
                    taxids[accession(line[1:].split("|", 1)[0])] = m.group(1)
    return taxids


class Species:
    """Maps a taxid to the taxid of its species, or None if it lies above any species."""

    def __init__(self, nodes_dmp):
        self.parent, self.rank, self.cache = {}, {}, {}
        with open(nodes_dmp, encoding="ascii") as f:
            for line in f:
                cols = line.split("\t|\t", 3)
                self.parent[cols[0]] = cols[1]
                self.rank[cols[0]] = cols[2]

    def of(self, taxid):
        if taxid not in self.cache:
            path, t = [], taxid
            result = None
            while t in self.parent and t not in self.cache:
                path.append(t)
                if self.rank[t] == "species":
                    result = t
                    break
                if self.parent[t] == t:
                    break
                t = self.parent[t]
            else:
                result = self.cache.get(t)
            for p in path:
                self.cache[p] = result
            self.cache.setdefault(taxid, result)
        return self.cache[taxid]


def read_rows(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=";"))


def read_tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def normalise(values):
    total = sum(values.values())
    return {k: v / total for k, v in values.items()} if total > 0 else dict(values)


def l1(estimate, truth, non_genome_mass=0.0):
    keys = set(estimate) | set(truth)
    return sum(abs(estimate.get(k, 0) - truth.get(k, 0)) for k in keys) + non_genome_mass


def topic_sets_of(specs):
    sets = []
    for spec in specs:
        label, _, path = spec.partition("=")
        if not path:
            sys.exit("--topics wants LABEL=CSV, got " + spec)
        sets.append((label, read_rows(path)))
    return sets


def match_total(match):
    return int(next((r["reads"] for r in match if r["name"] == "TOTAL"), 0) or 0)


def genome_level(args, reads, match, topic_sets):
    """Truth and estimates on genome nodes. Returns truth (read, org), the methods and info lines."""
    lengths = fasta_lengths(args.genomes)
    missing = [g for g in reads if g not in lengths]
    if missing:
        sys.exit("No length for truth genomes " + ", ".join(missing) + "; pass their fasta via --genomes.")
    truth_read = normalise(reads)
    truth_org = normalise({g: n / lengths[g] for g, n in reads.items()})
    methods, info = naive_genome_level("naive", "Genestrip", match, lengths, topic_sets)
    if args.ftmatch:
        ft_methods, ft_info = naive_genome_level("FT naive", "Genestrip-FT", read_rows(args.ftmatch), lengths,
                                                 topic_sets)
        methods += ft_methods
        info += ft_info
    for label, rows in topic_sets:
        tm_read, tm_org, off_read, off_org = {}, {}, 0.0, 0.0
        for row in rows:
            share, abundance = float(row["kmer share"]), float(row["abundance"])
            if row["rank"] == "GENOME":
                tm_read[row["name"]] = share
                tm_org[row["name"]] = abundance
            else:  # mass at a taxon topic: no genome gets it
                off_read += share
                off_org += abundance
        methods.append((label, tm_read, tm_org, off_read, off_org))
    return truth_read, truth_org, methods, info


def naive_genome_level(label, db, match, lengths, topic_sets):
    """The two naive variants (renormalised and strict) from one match CSV. Returns the methods and info lines."""
    genome_reads = {}
    other_reads = 0
    for row in match:
        if row["name"] == "TOTAL" or not row["reads"]:
            continue
        n = int(row["reads"])
        if row["rank"] == "GENOME":
            genome_reads[row["name"]] = genome_reads.get(row["name"], 0) + n
        else:
            other_reads += n
    genome_total = sum(genome_reads.values())
    classified = genome_total + other_reads
    dropped = other_reads / classified if classified else 0.0  # classified reads not at a genome node

    # Genome lengths for the naive organism share: fasta, else the topic model's expected hold at the genome
    # node (its k-mer positions, i.e. length - k + 1), else the genome's reads count as error.
    approx_len = dict(lengths)
    for _, rows in topic_sets:
        for row in rows:
            if row["rank"] == "GENOME" and row["name"] not in approx_len:
                approx_len[row["name"]] = float(row["expected hold"])
    unknown_len = {g: n for g, n in genome_reads.items() if g not in approx_len}
    unknown = sum(unknown_len.values()) / genome_total if genome_total else 0.0  # share of genome reads
    known_org = normalise({g: n / approx_len[g] for g, n in genome_reads.items() if g in approx_len})

    def scaled(shares, factor):
        return {g: v * factor for g, v in shares.items()}

    # Renormalised: the genome reads sum to one; only the unknown-length reads are lost (organism share).
    renorm_read = normalise(genome_reads)
    renorm_org = scaled(known_org, 1 - unknown)
    # Strict: shares of all classified reads; the dropped reads are mass not at a genome, as for the topic model.
    strict_read = {g: n / classified for g, n in genome_reads.items()} if classified else {}
    strict_org = scaled(known_org, (1 - dropped) * (1 - unknown))
    methods = [(label + " renorm.", renorm_read, renorm_org, 0.0, unknown),
               (label, strict_read, strict_org, dropped, dropped + (1 - dropped) * unknown)]

    total_reads = match_total(match)
    info = [f"{db}: reads {total_reads}, classified: {classified} "
            f"({classified / total_reads if total_reads else 0:.2%}), "
            f"of these not at a genome node: {other_reads} ({dropped:.2%})"]
    if unknown_len:
        info.append(f"{label} organism share: {len(unknown_len)} genomes without known length "
                    f"({sum(unknown_len.values())} reads) count as error.")
    return methods, info


def species_level(args, reads, match, topic_sets):
    """Truth and estimates summed up to species. Returns the truth read share, the methods (without organism
    share), info lines and species names where known."""
    if not args.taxonomy:
        sys.exit("--level species needs --taxonomy (nodes.dmp).")
    if not args.library and not all(g.isdigit() for g in reads):
        sys.exit("--level species needs --library (Kraken2 library fasta) to map accessions to tax ids.")
    species = Species(args.taxonomy)
    lib = library_taxids(args.library) if args.library else {}
    names = {}

    def of_genome(g):
        if g.isdigit():  # the genome is named after its tax id, see genome_of
            return species.of(g)
        if g not in lib:
            sys.exit(f"Genome {g} is not in the library {args.library}, so its species is unknown.")
        return species.of(lib[g])

    truth = Counter()
    for g, n in reads.items():
        s = of_genome(g)
        if s is None:
            sys.exit(f"Truth genome {g} (taxid {lib[g]}) has no species in the taxonomy.")
        truth[s] += n
    truth_read = normalise(truth)

    def two_variants(label, at_species, classified):
        """`renorm.` and strict shares from reads per species and all classified reads."""
        usable = sum(at_species.values())
        off = 1 - usable / classified if classified else 0.0
        strict = {s: n / classified for s, n in at_species.items()} if classified else {}
        return [(label + " renorm.", normalise(at_species), {}, 0.0, None),
                (label, strict, {}, off, None)], off

    def naive(label, db, rows):
        """Genestrip's reads at or below a species; a genome node's taxon is its parent taxid. REFINED nodes
        of Genestrip-FT carry an internal id, not a tax id, and count as above species."""
        at_species = Counter()
        classified = 0
        for row in rows:
            if row["name"] == "TOTAL" or not row["reads"]:
                continue
            n = int(row["reads"])
            classified += n
            if row["rank"] == "REFINED":
                continue
            s = species.of(row["parent taxid"] if row["rank"] == "GENOME" else row["taxid"])
            if s is not None:
                at_species[s] += n
            if row["rank"] == "species":
                names[row["taxid"]] = row["name"]
        found, dropped = two_variants(label, at_species, classified)
        total_reads = match_total(rows)
        return found, [f"{db}: reads {total_reads}, classified {classified} "
                       f"({classified / total_reads if total_reads else 0:.2%}), above species {dropped:.2%}"]

    methods, info = naive("naive", "Genestrip", match)
    if args.ftmatch:
        ft_methods, ft_info = naive("FT naive", "Genestrip-FT", read_rows(args.ftmatch))
        methods += ft_methods
        info += ft_info

    for label, rows in topic_sets:
        tm_read, off_read = Counter(), 0.0
        for row in rows:
            share = float(row["kmer share"])
            if row["rank"] == "REFINED":  # internal id, not a tax id: above species
                s = None
            elif row["rank"] == "GENOME":
                s = of_genome(row["name"])
            else:
                s = species.of(row["taxid"].lstrip("0"))
            if s is None:
                off_read += share
            else:
                tm_read[s] += share
        methods.append((label, dict(tm_read), {}, off_read, None))

    kraken_classified = None
    if args.kreport:
        cols = ["percent", "clade", "direct", "rank", "taxid", "name"]
        at_species = Counter()
        kraken_classified = 0
        unclassified = 0
        with open(args.kreport, encoding="utf-8") as f:
            for line in f:
                row = dict(zip(cols, line.rstrip("\n").split("\t")))
                name = row["name"].strip()
                if row["taxid"] == "1":
                    kraken_classified = int(row["clade"])
                elif row["taxid"] == "0":
                    unclassified = int(row["clade"])
                if row["rank"] == "S":
                    at_species[row["taxid"]] += int(row["clade"])
                    names[row["taxid"]] = name
        kraken, dropped = two_variants("Kraken2", at_species, kraken_classified)
        methods += kraken
        total = kraken_classified + unclassified
        info.append(f"Kraken2: read pairs {total}, classified {kraken_classified} "
                    f"({kraken_classified / total if total else 0:.2%}), above species {dropped:.2%}")

    if args.bracken:
        if kraken_classified is None:
            sys.exit("--bracken needs --kreport (the Kraken2 report Bracken ran on) for the classified reads.")
        rows = read_tsv(args.bracken)
        if any(r["taxonomy_lvl"] != "S" for r in rows):
            sys.exit("--bracken wants a species level Bracken output (bracken -l S).")
        est = {r["taxonomy_id"]: float(r["new_est_reads"]) for r in rows}
        for r in rows:
            names[r["taxonomy_id"]] = r["name"]
        own = {r["taxonomy_id"]: float(r["fraction_total_reads"]) for r in rows}
        usable = sum(est.values())
        off = 1 - usable / kraken_classified if kraken_classified else 0.0
        methods.append(("Bracken renorm.", own, {}, 0.0, None))
        methods.append(("Bracken", {s: n / kraken_classified for s, n in est.items()}, {}, off, None))
        info.append(f"Bracken: {usable:.0f} read pairs at species, not redistributed {off:.2%}")
    return truth_read, methods, info, names


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--level", choices=["genome", "species"], default="genome", help="level to compare on")
    ap.add_argument("--fastq", required=True, help="ISS fastq (R1) of the sample, for the ground truth")
    ap.add_argument("--genomes", nargs="+", help="fasta file(s) of the simulated genomes, for lengths "
                                                 "(genome level)")
    ap.add_argument("--match", required=True, help="Genestrip match CSV of the sample")
    ap.add_argument("--ftmatch", help="Genestrip-FT ftmatch CSV of the sample; adds the naive FT variants")
    ap.add_argument("--topics", nargs="*", default=[], metavar="LABEL=CSV", help="topicmodelcsv file(s) of the sample")
    ap.add_argument("--taxonomy", help="nodes.dmp of the database taxonomy (species level)")
    ap.add_argument("--library", help="Kraken2 library fasta with kraken:taxid headers, maps the genomes to "
                                      "taxids (species level)")
    ap.add_argument("--kreport", help="Kraken2 report of the sample (species level)")
    ap.add_argument("--bracken", help="Bracken output of the sample, run with -l S (species level)")
    ap.add_argument("--top", type=int, default=0, help="table rows: truth genomes plus the N largest false positives")
    ap.add_argument("--out", help="also write the per-genome table to this CSV")
    ap.add_argument("--sample", help="name of the sample in the summary and plot; default: from the fastq name")
    ap.add_argument("--summary-out", help="also write the per-method summary (L1 etc.) to this CSV")
    ap.add_argument("--plot", help="also draw estimate against truth per genome into this file (.pdf, .png)")
    args = ap.parse_args()
    sample = args.sample or re.sub(r"(_R?1)?\.f(ast)?q(\.gz)?$", "", args.fastq.replace("\\", "/").split("/")[-1])
    if args.level == "genome" and (args.kreport or args.bracken):
        sys.exit("Kraken2 and Bracken have no genome nodes; compare them with --level species.")
    if args.level == "genome" and not args.genomes:
        sys.exit("--level genome needs --genomes for the genome lengths.")

    reads = truth_reads(args.fastq)
    if not reads:
        sys.exit("No reads in " + args.fastq + ", so there is no ground truth to compare against.")
    match = read_rows(args.match)
    topic_sets = topic_sets_of(args.topics)
    names = {}
    if args.level == "genome":
        truth_read, truth_org, methods, info = genome_level(args, reads, match, topic_sets)
    else:
        truth_read, methods, info, names = species_level(args, reads, match, topic_sets)
        truth_org = None
    unit = args.level

    # Report.
    print(f"Sample: {args.fastq} ({unit} level)")
    for line in info:
        print(line)
    print()
    print(f"| method | L1 read share | L1 organism share | mass not at a {unit} (read / organism) "
          "| FP mass (read / organism) |")
    print("|---|---|---|---|---|")
    summary = []

    def fmt(v):
        return "-" if v is None else f"{v:.4f}"

    for label, est_read, est_org, off_read, off_org in methods:
        fp_read = sum(v for g, v in est_read.items() if g not in truth_read)
        l1_read = l1(est_read, truth_read, off_read)
        if truth_org is None:
            fp_org = l1_org = None
        else:
            fp_org = sum(v for g, v in est_org.items() if g not in truth_org)
            l1_org = l1(est_org, truth_org, off_org)
        summary.append([sample, label, l1_read, "" if l1_org is None else l1_org, off_read,
                        "" if off_org is None else off_org, fp_read, "" if fp_org is None else fp_org])
        print(f"| {label} | {l1_read:.4f} | {fmt(l1_org)} "
              f"| {off_read:.4f} / {fmt(off_org)} | {fp_read:.4f} / {fmt(fp_org)} |")
    print()
    if args.summary_out:
        with open(args.summary_out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter=";")
            w.writerow(SUMMARY_HEADER)
            w.writerows(summary)

    false_pos = sorted({g for _, r, _, _, _ in methods for g in r if g not in truth_read},
                       key=lambda g: -max(m[1].get(g, 0) for m in methods))
    shown = sorted(truth_read, key=lambda g: -truth_read[g]) + false_pos[:args.top]
    if truth_org is None:
        header = ["species", "name", "truth read"]
        for label, *_ in methods:
            header.append(label + " read")
        table = [[g, names.get(g, "")] + [truth_read.get(g, 0)] + [m[1].get(g, 0) for m in methods]
                 for g in shown]
        text_cols = 2
    else:
        header = ["genome", "truth read", "truth org"]
        for label, *_ in methods:
            header += [label + " read", label + " org"]
        table = []
        for g in shown:
            row = [g, truth_read.get(g, 0), truth_org.get(g, 0)]
            for _, est_read, est_org, _, _ in methods:
                row += [est_read.get(g, 0), est_org.get(g, 0)]
            table.append(row)
        text_cols = 1
    print("| " + " | ".join(header) + " |")
    print("|" + "---|" * len(header))
    for row in table:
        print("| " + " | ".join(row[:text_cols]) + " | " + " | ".join(f"{v:.4f}" for v in row[text_cols:]) + " |")

    if args.out:
        with open(args.out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter=";")
            w.writerow(header)
            w.writerows(table)

    if args.plot:
        plot_scatter(args.plot, sample, truth_read, truth_org, methods, unit)


SUMMARY_HEADER = ["sample", "method", "l1 read", "l1 org", "not at genome read", "not at genome org",
                  "fp read", "fp org"]


def plot_scatter(path, sample, truth_read, truth_org, methods, unit="genome"):
    """Estimate against truth per genome (species): one column per method, one row per kind of share (small
    multiples, since the methods would hide each other in one panel), log-log so that rare genomes show. A
    point on the diagonal is exact. Values below a floor a decade under the rarest truth genome - zeros
    included, which a log axis cannot show - are drawn on the floor line: a truth genome a method missed lies
    on the bottom line, a false positive (truth 0) on the left one. Without organism shares (species level)
    only the read share row is drawn."""
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # plot_style.py sits beside this script
    import plot_style
    plt = plot_style.setup()
    rows = [(truth_read, 1, "Read-Anteil")]
    if truth_org is not None:
        rows.append((truth_org, 2, "Organismen-Anteil"))
    cols = min(len(methods), 4) if truth_org is None else len(methods)
    grid_rows = len(rows) * -(-len(methods) // cols)
    fig, axes = plt.subplots(grid_rows, cols, figsize=(2.5 * cols + 0.6, 2.8 * grid_rows), squeeze=False)
    for r, (truth, which, kind) in enumerate(rows):
        floor = min(v for v in truth.values() if v > 0) / 10
        top = 1.5
        for i, m in enumerate(methods):
            gr, gc = (r, i) if truth_org is not None else divmod(i, cols)
            ax = axes[gr][gc]
            est = m[which]
            genomes = set(truth) | {g for g, v in est.items() if v > 0}
            xs = [max(truth.get(g, 0), floor) for g in genomes]
            ys = [max(est.get(g, 0), floor) for g in genomes]
            ax.plot([floor, top], [floor, top], color=plot_style.MUTED, linewidth=1, linestyle="--", zorder=2)
            ax.axhline(floor, color=plot_style.GRID, linewidth=2, zorder=1)
            ax.axvline(floor, color=plot_style.GRID, linewidth=2, zorder=1)
            ax.scatter(xs, ys, s=30, edgecolors="white", linewidths=0.8, zorder=3, **plot_style.style(i))
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(floor / 1.6, top)
            ax.set_ylim(floor / 1.6, top)
            ax.set_aspect("equal")
            ax.set_title(m[0])
            if gc == 0:
                ax.set_ylabel(f"{kind}\ngeschätzt")
            else:
                ax.tick_params(labelleft=False)
            if gr == grid_rows - 1 or (truth_org is None and i + cols >= len(methods)):
                ax.set_xlabel("wahr")  # only under the lowest panel of a column, else it hits the title below
    for gr in range(grid_rows):  # empty cells of the last row
        for gc in range(cols):
            if not axes[gr][gc].has_data():
                axes[gr][gc].set_visible(False)
    what = "Genom" if unit == "genome" else "Spezies"
    fig.suptitle(f"Probe {sample}: Schätzung gegen Wahrheit pro {what}. Diagonale = exakt; "
                 f"Punkte auf der unteren/linken Linie: 0 oder fast 0", color=plot_style.INK, fontsize=9.5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print(f"Plot written to {path}")


if __name__ == "__main__":
    main()
