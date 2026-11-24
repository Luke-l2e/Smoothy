#!/usr/bin/env python3
"""Builds the sample overview page docs/proben-taxtree/index.html.

Reads the genome path counts of the protozoa DB (normal and refined "ft" tree) and the
genome mix of the InSilicoSeq run, and builds one JSON with both trees, every genome's
k-mer path and the planned test samples. The JSON is embedded into the template
sample_tree_page.html, so the resulting index.html opens as a plain file in any browser.
Edit the planned samples in SAMPLES below.

usage: python3 sample_tree_page.py DATA_DIR OUT_DIR      (e.g. data docs/proben-taxtree)
"""
import csv, json, os, sys
from collections import defaultdict

DATA, OUT = sys.argv[1], sys.argv[2]
PROJ = f"{DATA}/projects/protozoa"

# Planned samples, by genome taxid. Proportions: None = equal shares.
SAMPLES = [
    {"id": "p1", "label": "1 · Weit entfernt", "title": "Probe 1: weit entfernte Gattungen (Kontrolle)",
     "text": "Je ein Genom aus Gattungen ohne nahe Verwandte in der DB. Fast alle k-mers liegen am Genomknoten, "
             "also sollten alle Methoden fast perfekt sein. Klappt das nicht, liegt der Fehler im Ablauf.",
     "taxids": ["5741", "484906", "5763", "944170", "508771", "310047"]},
    {"id": "p2", "label": "2 · Eng verwandt, alle drin", "title": "Probe 2: eng verwandt, alle Verwandten in der Probe",
     "text": "Paare und Gruppen, deren k-mers großteils an einem gemeinsamen Knoten liegen, jeweils komplett in der Probe. "
             "Die Methoden müssen die geteilten k-mers auf Genome aufteilen, die wirklich in der Probe sind.",
     "taxids": ["5858", "5824", "5661", "435258", "185431", "679716", "5694"]},
    {"id": "p3a", "label": "3a · Nur ein Verwandter", "title": "Probe 3a: k-mers weit oben, nur ein Verwandter in der Probe",
     "text": "Aus jeder Gruppe von Probe 2 nur ein Genom. Die k-mers liegen hoch im Baum, weil der Verwandte in der DB ist, "
             "nicht weil er in der Probe ist. Naive Methoden verlieren die Reads an den höheren Knoten.",
     "taxids": ["435258", "5824", "185431"]},
    {"id": "p3b", "label": "3b · Gegenstück", "title": "Probe 3b: Gegenstück zu 3a",
     "text": "Jeweils der andere Partner aus Probe 2. Zusammen mit 3a ist jedes Genom einmal allein in der Probe.",
     "taxids": ["5661", "5858", "679716"]},
    {"id": "p5", "label": "5 · Leave-one-out", "title": "Probe 5: Genom fehlt in der DB (leave-one-out)",
     "text": "Je ein eigenes Projekt (protozoa_loo_<taxid>), dessen DB genau eines dieser Genome auslässt; die Reads "
             "stammen trotzdem von ihm. L. infantum (Verwandter L. donovani bleibt in der DB), C. parvum (C. hominis "
             "bleibt) und Babesia bovis (kein naher Verwandter). Der Baum zeigt die volle protozoa-DB: die markierten "
             "Pfade sind die, die in der jeweiligen LOO-DB fehlen bzw. anders aussehen.",
     "taxids": ["435258", "353152", "484906"]},
    {"id": "mix", "label": "Realistische Mischung", "title": "Realistische Mischung (bisheriger ISS-Mix)",
     "text": "Der bisherige Mix aus insilicoseq/mix.tsv, mit Anteilen über etwa zwei Größenordnungen.",
     "taxids": None},
]
FT_ONLY_NOTE = "REFINED"


def read_mix():
    out = []
    for line in open(f"{PROJ}/insilicoseq/mix.tsv", encoding="utf-8"):
        line = line.split("#")[0].strip()
        if line:
            t, w = line.split()[:2]
            out.append((t, float(w)))
    tot = sum(w for _, w in out)
    return [(t, w / tot) for t, w in out]


def read_tree(path):
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter=";"))[1:]
    by_genome = defaultdict(list)
    for r in rows:
        if ":" in r[0]:          # k-mers of taxa outside the project, not genomes of the DB
            continue
        by_genome[r[0]].append(r)
    nodes, genomes = {}, {}
    for key, rs in by_genome.items():
        rs.sort(key=lambda r: int(r[2]))
        total = sum(int(r[6]) for r in rs)
        path = []
        for r in rs:
            nid = "g" + key if r[5] == "GENOME" else r[3]
            if nid not in nodes:
                nodes[nid] = {"name": r[4], "rank": r[5], "parent": None}
            path.append([nid, int(r[6]), float(r[7])])
        for (child, *_), (par, *_) in zip(path, path[1:]):
            old = nodes[child]["parent"]
            if old not in (None, par):
                print(f"warning: {path} {child} has parents {old} and {par}", file=sys.stderr)
            nodes[child]["parent"] = par
        taxon = rs[1]
        genomes[key] = {"taxid": rs[0][1], "name": taxon[4], "rank": taxon[5], "kmers": total, "path": path}
        nodes["g" + key]["name"] = taxon[4]
    return nodes, genomes


std_nodes, std_genomes = read_tree(f"{PROJ}/db/protozoa_kmergenomepathcountscsv.csv")
ft_nodes, ft_genomes = read_tree(f"{PROJ}/db/protozoa_ftkmergenomepathcountscsv.csv")
key_of = {g["taxid"]: k for k, g in std_genomes.items()}
samples = []
for s in SAMPLES:
    s = dict(s)
    if s["taxids"] is None:
        mix = read_mix()
        s["genomes"] = [key_of[t] for t, _ in mix]
        s["props"] = [round(w, 4) for _, w in mix]
    else:
        s["genomes"] = [key_of[t] for t in s["taxids"]]
        s["props"] = None
    del s["taxids"]
    samples.append(s)
data = {"trees": {"std": {"nodes": std_nodes, "genomes": std_genomes},
                  "ft": {"nodes": ft_nodes, "genomes": ft_genomes}},
        "samples": samples}
blob = json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
template = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "sample_tree_page.html"), encoding="utf-8").read()
os.makedirs(OUT, exist_ok=True)
with open(f"{OUT}/index.html", "w", encoding="utf-8") as f:
    f.write(template.replace("__DATA__", blob, 1))
print("written", f"{OUT}/index.html", len(std_genomes), "genomes,", len(std_nodes), "/", len(ft_nodes), "nodes")
