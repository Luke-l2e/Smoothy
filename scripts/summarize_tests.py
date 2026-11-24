#!/usr/bin/env python3
"""Summarises the runs of run_tests.sh over their seeds.

Reads <OUT>/runs/*/{genome,species}_summary.csv (one row per method, written by
compare_abundances.py) and the run's `done` file (sample file, error level, seed), and writes per
level <OUT>/summary_<level>.csv with one row per sample file, error level and method: number of
seeds, and mean and standard deviation of every error column. Also collects the run times of the
steps (runs/*/timings.csv) into timings_runs.csv and timings_summary.csv.

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

# Run times per step (runs/*/timings.csv, written by run_tests.sh): one row per run and step, and
# mean and standard deviation per error level and step. smoothy's times are per batch and stay in
# smoothy_timings.csv.
rows, by_step = [], defaultdict(list)
for path in sorted(glob.glob(f"{OUT}/runs/*/timings.csv")):
    run = os.path.dirname(path)
    if not os.path.exists(f"{run}/done"):
        continue
    meta = dict(line.strip().split("=", 1) for line in open(f"{run}/done") if "=" in line)
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter=";"):
            rows.append([os.path.basename(run), meta["sample"], meta["error"], meta["seed"],
                         meta.get("n_reads", ""), row["step"], row["seconds"]])
            by_step[(meta["error"], meta.get("n_reads", ""), row["step"])].append(float(row["seconds"]))
if rows:
    with open(f"{OUT}/timings_runs.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["run", "sample", "error", "seed", "n_reads", "step", "seconds"])
        w.writerows(rows)
    with open(f"{OUT}/timings_summary.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["error", "n_reads", "step", "runs", "seconds mean", "seconds sd"])
        for key in sorted(by_step):
            v = by_step[key]
            w.writerow(list(key) + [len(v), f"{statistics.mean(v):.1f}",
                                    f"{statistics.stdev(v):.1f}" if len(v) > 1 else ""])
    print(f"written {OUT}/timings_runs.csv and timings_summary.csv: {len(rows)} rows")
