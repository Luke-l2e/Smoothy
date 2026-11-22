#!/usr/bin/env python3
"""Prepares an InSilicoSeq run over a chosen set of genomes from a Genestrip project.

The genomes come from the per-accession FASTA files that the goal `extractrefseqcsv` writes into
data/projects/<db>/fasta, so every simulated read has its counterpart in the database. Their headers
carry the tax id (`>ACC|kraken:taxid|TAXID`), and InSilicoSeq keeps the record id in the read name,
so the ground truth can be read off every read.

InSilicoSeq gives each FASTA record its own abundance, not each genome, and protozoan genomes consist
of up to thousands of contigs. The abundance file written here therefore splits each genome's read
share over its records in proportion to their length, so that the genome as a whole gets the
requested share of the reads. Records shorter than `--min-length` are left out (InSilicoSeq skips
records shorter than its fragment anyway), and the share of the remaining ones is rescaled.

Input is a tab-separated file with one genome per line: tax id and read share (any positive numbers,
they are normalised). Lines starting with '#' are ignored.

Output (prefix given by `--out`):
    <out>_genomes.fna     the records to simulate from, for `iss generate --genomes`
    <out>_abundance.txt   record id and abundance, for `iss generate --abundance_file`
    <out>_truth.tsv       per genome: tax id, records, length, requested read share

Example:
    python iss_sample.py --fasta ../data/projects/protozoa/fasta --genomes protozoa_mix.tsv \
        --out ../data/projects/protozoa/insilicoseq/mix
"""
import argparse
import os
import sys
from collections import defaultdict


def read_selection(path):
    shares = {}
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                taxid, share = line.split()[:2]
                shares[taxid] = float(share)
    total = sum(shares.values())
    return {t: s / total for t, s in shares.items()}


def read_record(path):
    """Returns (record id, tax id, sequence lines) of a single-record FASTA file."""
    with open(path) as f:
        header = f.readline()
        lines = f.readlines()
    rid = header[1:].split(None, 1)[0]
    return rid, rid.rsplit("|", 1)[-1], header, lines


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="directory with the files extractrefseqcsv wrote")
    parser.add_argument("--genomes", required=True, help="tab-separated tax id and read share per genome")
    parser.add_argument("--out", required=True, help="prefix of the output files")
    parser.add_argument("--min-length", type=int, default=1000, help="shortest record to keep (default 1000)")
    args = parser.parse_args()

    shares = read_selection(args.genomes)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)

    lengths = defaultdict(dict)  # taxid -> record id -> length
    short = defaultdict(int)
    with open(args.out + "_genomes.fna", "w") as out:
        for name in sorted(os.listdir(args.fasta)):
            path = os.path.join(args.fasta, name)
            with open(path) as f:
                header = f.readline()
            rid = header[1:].split(None, 1)[0]
            taxid = rid.rsplit("|", 1)[-1]
            if taxid not in shares:
                continue
            rid, taxid, header, lines = read_record(path)
            length = sum(len(line.strip()) for line in lines)
            if length < args.min_length:
                short[taxid] += 1
                continue
            lengths[taxid][rid] = length
            out.write(header)
            out.writelines(lines)

    missing = [t for t in shares if t not in lengths]
    if missing:
        sys.exit(f"no records of at least {args.min_length} bp for tax id(s) {', '.join(missing)}")

    with open(args.out + "_abundance.txt", "w") as ab, open(args.out + "_truth.tsv", "w") as truth:
        truth.write("taxid\trecords\tlength\tshort_records_dropped\tread_share\n")
        for taxid, share in shares.items():
            total = sum(lengths[taxid].values())
            for rid, length in lengths[taxid].items():
                ab.write(f"{rid}\t{share * length / total:.10g}\n")
            truth.write(f"{taxid}\t{len(lengths[taxid])}\t{total}\t{short[taxid]}\t{share:.6f}\n")
            print(f"{taxid}: {len(lengths[taxid])} records, {total / 1e6:.1f} Mbp, share {share:.4f}"
                  + (f", {short[taxid]} short records dropped" if short[taxid] else ""), file=sys.stderr)
            if total < 1_000_000:
                print(f"  WARNING: {taxid} has only {total} bp - likely an organelle or single gene, not a genome",
                      file=sys.stderr)


if __name__ == "__main__":
    main()
