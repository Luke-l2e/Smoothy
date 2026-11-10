#!/usr/bin/env python3
"""Minimal version of compare_abundances.py: genome level, read shares only, no plots.

For one simulated sample it compares these read-share vectors against the ground truth:

* truth  - reads per genome, counted from the InSilicoSeq read names (@NC_040360.1_0_0/1).
* naive  - Genestrip's classified reads at genome nodes, divided by all classified reads.
           Reads Genestrip places higher in the tree (genus, family ...) reach no genome -> error.
* naive renorm. - the same reads at genome nodes, divided by the reads at genome nodes only: the reads
           above GENOME are ignored (in effect spread proportionally over the genomes).
* topic  - `kmer share` (theta) of the genome topics from a topicmodelcsv file.
           Mass at taxon topics reaches no genome -> error.

Distance: L1 = sum_g |estimate_g - truth_g| + mass not at a genome   (0 = perfect, 2 = disjoint).

Example:
    python compare.py ../data/projects/viral/fastq/miseq_cov_R1.fastq \
        ../data/projects/viral/eval/weighted/viral_match_miseq_cov.csv \
        ../data/projects/viral/eval/weighted/viral_topicmodelcsv_miseq_cov.csv
"""
import csv
import re
import sys
from collections import Counter


def accession(seq_id):
    return seq_id.split(".")[0]  # NC_040360.1 -> NC_040360


def truth_shares(fastq):
    counts = Counter()
    with open(fastq, encoding="ascii") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:  # header line of each fastq record
                source = re.match(r"^@(.+)_\d+_\d+(/[12])?$", line.strip()).group(1)
                counts[accession(source)] += 1
    total = sum(counts.values())
    return {g: n / total for g, n in counts.items()}


def read_csv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=";"))
    return None


def naive_shares(match_csv):
    """Returns (shares per genome, share of classified reads not at a genome)."""
    at_genome, classified = Counter(), 0
    for row in read_csv(match_csv):
        if row["name"] == "TOTAL" or not row["reads"]:
            continue
        classified += int(row["reads"])
        if row["rank"] == "GENOME":
            at_genome[row["name"]] += int(row["reads"])
    shares = {g: n / classified for g, n in at_genome.items()}
    return shares, 1 - sum(shares.values())


def naive_renorm_shares(match_csv):
    """Like naive_shares, but divided by the reads at genomes only: reads above GENOME are ignored."""
    shares, off = naive_shares(match_csv)
    return {g: v / (1 - off) for g, v in shares.items()}, 0.0


def topic_shares(topic_csv):
    """Returns (theta per genome, theta mass at taxon topics)."""
    shares, off = {}, 0.0
    for row in read_csv(topic_csv):
        if row["rank"] == "GENOME":
            shares[row["name"]] = float(row["kmer share"])
        else:
            off += float(row["kmer share"])
    return shares, off


def l1(estimate, truth, off):
    genomes = set(estimate) | set(truth)
    return sum(abs(estimate.get(g, 0) - truth.get(g, 0)) for g in genomes) + off


def main():
    if len(sys.argv) != 4:
        sys.exit("usage: compare.py SAMPLE_R1.fastq MATCH.csv TOPICMODEL.csv")
    truth = truth_shares(sys.argv[1])
    methods = {"naive": naive_shares(sys.argv[2]),
               "naive renorm.": naive_renorm_shares(sys.argv[2]),
               "topic": topic_shares(sys.argv[3])}

    print("| method | L1 | mass not at a genome |")
    print("|---|---|---|")
    for name, (est, off) in methods.items():
        print(f"| {name} | {l1(est, truth, off):.4f} | {off:.4f} |")
    print()
    print("| genome | truth | " + " | ".join(methods) + " |")
    print("|---|---|" + "---|" * len(methods))
    for g in sorted(truth, key=truth.get, reverse=True):
        print(f"| {g} | {truth[g]:.4f} | " + " | ".join(f"{est.get(g, 0):.4f}" for est, _ in methods.values()) + " |")


if __name__ == "__main__":
    main()
