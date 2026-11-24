#!/usr/bin/env python3
"""Builds the topic result page docs/topic-ergebnisse/index.html.

For every test run of run_tests.sh in the given projects it collects the four topic model fits (standard
and refined "ft" tree, with and without consistency weight, i.e. weighted/ and unweighted/), the ground
truth per genome (genome_comparison.csv) and the per-method summary (genome_summary.csv), together with
both trees of the project's database (k-mer genome path counts). Everything is embedded as one JSON into
the template topic_result_page.html, so index.html opens as a plain file in any browser.

The page shows, per run and fit:
  * a bar chart of the topics, sorted by share or in taxonomic order, with the truth marked, and
  * a tree that has only the topics as nodes (each topic hangs below its nearest ancestor that is a topic
    too), coloured from cold (blue) to hot (red) by its share.

usage: python3 topic_result_page.py DATA_DIR OUT_DIR PROJECT...
       (e.g. data docs/topic-ergebnisse protozoa protozoa_loo_435258 protozoa_loo_353152)
"""
import csv
import glob
import json
import os
import re
import sys

DATA, OUT, PROJECTS = sys.argv[1], sys.argv[2], sys.argv[3:]
FITS = {"std-w": ("weighted", "topicmodelcsv", "weighted"),
        "std-u": ("unweighted", "topicmodelcsv", "unweighted"),
        "ft-w": ("weighted", "fttopicmodelcsv", "FT weighted"),
        "ft-u": ("unweighted", "fttopicmodelcsv", "FT unweighted")}
RUN = re.compile(r"^(?P<sample>.+)_(?P<error>perfect|err\d+)_s(?P<seed>\d+)$")


def rows(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.reader(f, delimiter=";"))


def read_tree(path):
    """Node id -> [name, rank, parent id]; also genome key -> (genome tax id, organism name) for stored genomes.

    Every genome's row lists its path from the genome (position 0) up to the root, so each node's parent is the
    next node of any path through it. Genomes of the RefSeq LCA update that the DB does not store ("2759:100879")
    are walked too, since topics can sit at nodes only they reach."""
    paths = {}
    for r in rows(path)[1:]:
        paths.setdefault(r[0], []).append(r)
    nodes, genomes = {}, {}
    for key, rs in paths.items():
        rs.sort(key=lambda r: int(r[2]))
        for child, parent in zip(rs, rs[1:] + [None]):
            entry = nodes.setdefault(child[3], [child[4], child[5], None])
            if parent is not None:
                entry[2] = parent[3]
        if rs[0][8] == "true" and rs[0][5] == "GENOME":
            # the organism is the first taxon above the genome node; in the FT tree refined nodes may come first
            taxon = next((r for r in rs[1:] if r[5] != "REFINED"), None)
            genomes[rs[0][3]] = (rs[0][1], taxon[4] if taxon else rs[0][1])
    return nodes, genomes


def num(s):
    return float(s) if s not in ("", None) else None


def compact(x):
    return float(f"{x:.6g}")


def excluded_taxids(project_dir):
    path = f"{project_dir}/taxids.txt"
    if not os.path.exists(path):
        return []
    return [m.group(1) for line in open(path, encoding="utf-8") if (m := re.match(r"\s*-(\d+)", line))]


data = {"projects": {}, "taxidNames": {}}  # taxidNames: genome tax id -> name, to name leave-one-out genomes
for project in PROJECTS:
    pdir = f"{DATA}/projects/{project}"
    std_nodes, std_genomes = read_tree(f"{pdir}/db/{project}_kmergenomepathcountscsv.csv")
    ft_nodes, ft_genomes = read_tree(f"{pdir}/db/{project}_ftkmergenomepathcountscsv.csv")
    names = {k: v[1] for k, v in {**std_genomes, **ft_genomes}.items()}
    data["taxidNames"].update({v[0]: v[1] for v in {**std_genomes, **ft_genomes}.values()})
    key_of_taxid = {v[0]: k for k, v in {**std_genomes, **ft_genomes}.items()}
    runs = []
    for run_dir in sorted(glob.glob(f"{pdir}/tests/runs/*")):
        run = os.path.basename(run_dir)
        m = RUN.match(run)
        if not m or not os.path.exists(f"{run_dir}/done"):
            continue
        fits = {}
        for fit, (sub, goal, _) in FITS.items():
            path = f"{run_dir}/{sub}/{project}_{goal}_{run}.csv"
            if os.path.exists(path):
                fits[fit] = [[r[0], int(r[3]), int(r[4]), compact(float(r[5])), compact(float(r[6])),
                              compact(float(r[7])), compact(float(r[8]))] for r in rows(path)[1:] if r]
        if not fits:
            continue
        truth, missing = {}, {}
        cmp_path = f"{run_dir}/genome_comparison.csv"
        if os.path.exists(cmp_path):
            with open(cmp_path, encoding="utf-8") as f:
                for r in csv.DictReader(f, delimiter=";"):
                    t_read, t_org = num(r["truth read"]), num(r["truth org"])
                    if not t_read:
                        continue
                    key = key_of_taxid.get(r["genome"])
                    target = truth if key else missing
                    target[key or r["genome"]] = [compact(t_read), compact(t_org) if t_org is not None else None]
        summary = {}
        sum_path = f"{run_dir}/genome_summary.csv"
        if os.path.exists(sum_path):
            with open(sum_path, encoding="utf-8") as f:
                by_method = {r["method"]: r for r in csv.DictReader(f, delimiter=";")}
            for fit, (_, _, method) in FITS.items():
                r = by_method.get(method)
                if r:
                    summary[fit] = [compact(num(r[c])) if num(r[c]) is not None else None
                                    for c in ("l1 read", "l1 org", "not at genome read", "not at genome org",
                                              "fp read", "fp org")]
        runs.append({"id": run, **m.groupdict(), "truth": truth, "missing": missing, "summary": summary,
                     "fits": fits})
    if not runs:
        print(f"{project}: no finished runs, skipped", file=sys.stderr)
        continue
    data["projects"][project] = {"excluded": excluded_taxids(pdir), "names": names,
                                 "trees": {"std": std_nodes, "ft": ft_nodes}, "runs": runs}
    print(f"{project}: {len(runs)} runs, {len(std_nodes)} / {len(ft_nodes)} nodes")

blob = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "topic_result_page.html")
template = open(template_path, encoding="utf-8").read()
os.makedirs(OUT, exist_ok=True)
with open(f"{OUT}/index.html", "w", encoding="utf-8") as f:
    f.write(template.replace("__DATA__", blob, 1))
print(f"written {OUT}/index.html ({os.path.getsize(f'{OUT}/index.html') / 1e6:.1f} MB)")
