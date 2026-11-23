#!/usr/bin/env python3
"""Summarises the runs of run_tests.sh over their seeds.

Reads <OUT>/runs/*/{genome,species}_summary.csv (one row per method, written by
compare_abundances.py) and the run's `done` file (sample file, error level, seed), and writes per
level <OUT>/summary_<level>.csv with one row per sample file, error level and method: number of
seeds, and mean and standard deviation of every error column.

usage: python3 summarize_tests.py OUT_DIR
"""
import csv, glob, os, statistics, sys
from collections import defaultdict

OUT = sys.argv[1]
for level in ("genome", "species"):
    groups, columns, order = defaultdict(lambda: defaultdict(list)), None, []
    for path in sorted(glob.glob(f"{OUT}/runs/*/{level}_summary.csv")):
        run = os.path.dirname(path)
        if not os.path.exists(f"{run}/done"):
            continue
        meta = dict(line.strip().split("=", 1) for line in open(f"{run}/done") if "=" in line)
        with open(path, encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter=";"))
        for row in rows:
            columns = columns or [c for c in row if c not in ("sample", "method")]
            key = (meta["sample"], meta["error"], row["method"])
            if key not in groups:
                order.append(key)
            for c in columns:
                if row[c] != "":
                    groups[key][c].append(float(row[c]))
            groups[key]["_seeds"].append(meta["seed"])
    if not order:
        continue
    with open(f"{OUT}/summary_{level}.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["sample", "error", "method", "seeds"]
                   + [f"{c} {s}" for c in columns for s in ("mean", "sd")])
        for key in sorted(order):
            g = groups[key]
            out = list(key) + [len(g["_seeds"])]
            for c in columns:
                v = g.get(c, [])
                out += [f"{statistics.mean(v):.6g}" if v else "",
                        f"{statistics.stdev(v):.3g}" if len(v) > 1 else ""]
            w.writerow(out)
    print(f"written {OUT}/summary_{level}.csv: {len(order)} rows")
