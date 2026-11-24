#!/usr/bin/env python3
"""Draws the L1 distance to the truth for the test samples p1-p4, one figure per sample group.

Reads summary_genome.csv and summary_species.csv written by summarize_tests.py (mean and standard deviation
over the seeds) and writes <OUT>/plots/l1_<group>.png and .pdf. Each figure has three rows and two columns:

    rows    genome level, read share     naive / core topic model / FT topic model
            genome level, organism share naive / core topic model / FT topic model
            species level, read share    core topic model / FT topic model / Bracken
    columns perfect reads, err10 reads

Kraken2 and Bracken stop at the species, so they only appear on species level. The topic models are the
weighted ones (consistency weight on, the default), naive and Bracken the strict variants: mass a method
cannot give to a genome (species) counts as error, as it does for the topic models. Whiskers are +- one
standard deviation over the seeds (p4 has one seed, so none). The axis is logarithmic, since the values span
four orders of magnitude; values below FLOOR are drawn on the floor line.

usage: python3 plot_tests.py OUT_DIR SAMPLES_DIR     (e.g. data/projects/protozoa/tests
                                                      data/projects/protozoa/insilicoseq/samples)
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_style  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

OUT, SAMPLES = sys.argv[1], sys.argv[2]
FLOOR = 1e-4

# One colour and marker per method in every figure. Two of the three colour slots per panel are the
# palette's first slots; the naive baseline is neutral grey. Validated all-pairs for colour-vision deficiency
# (blue, orange, aqua); every method also has its own marker and a legend.
METHODS = {
    "naive": ("naive", "#8a8984", "s"),
    "weighted": ("Topic Model (core)", "#2a78d6", "o"),
    "FT weighted": ("Topic Model (FT)", "#eb6834", "D"),
    "Bracken": ("Bracken", "#1baf7a", "^"),
}
ROWS = [("genome", "l1 read", "Genomebene, Read-Anteil", ["naive", "weighted", "FT weighted"]),
        ("genome", "l1 org", "Genomebene, Organismen-Anteil", ["naive", "weighted", "FT weighted"]),
        ("species", "l1 read", "Artebene, Read-Anteil", ["weighted", "FT weighted", "Bracken"])]
ERRORS = [("perfect", "fehlerfreie Reads"), ("err10", "err10 (10 % Fehler pro Base)")]
GROUPS = {
    "p1": ("Probe 1: weit entfernte Gattungen (Kontrolle)", ["p1_equal", "p1_skew"]),
    "p2": ("Probe 2: eng verwandt, alle Verwandten in der Probe", ["p2_equal", "p2_skew_a", "p2_skew_b"]),
    "p3": ("Probe 3: nur ein Verwandter in der Probe", ["p3a_equal", "p3a_skew", "p3b_equal", "p3b_skew"]),
    "p4": ("Probe 4: je ein Genom allein (ein Seed)", None),
}


def read(level):
    with open(f"{OUT}/summary_{level}.csv", encoding="utf-8") as f:
        return {(r["sample"], r["error"], r["method"]): r for r in csv.DictReader(f, delimiter=";")}


def p4_label(sample):
    """p4_5741 -> 'Giardia duodenalis' from the comment in the sample file, else the tax id."""
    path = f"{SAMPLES}/{sample}.tsv"
    if os.path.exists(path):
        for line in open(path, encoding="utf-8"):
            if line.strip() and not line.startswith("#") and "#" in line:
                return line.split("#", 1)[1].strip()
    return sample.split("_", 1)[1]


data = {"genome": read("genome"), "species": read("species")}
plt = plot_style.setup()
os.makedirs(f"{OUT}/plots", exist_ok=True)

for group, (title, samples) in GROUPS.items():
    if samples is None:
        samples = sorted({s for s, _, _ in data["genome"] if s.startswith(group + "_")}, key=p4_label)
        labels = [p4_label(s) for s in samples]
    else:
        # p2_skew_a -> "skew a", p3a_equal -> "3a equal" (keeps the sub-sample when the group has several)
        labels = [(s.split("_", 1)[0][1:] + " " if s.split("_", 1)[0] != group else "")
                  + s.split("_", 1)[1].replace("_", " ") for s in samples]
    wide = len(samples) > 5
    fig, axes = plt.subplots(len(ROWS), len(ERRORS), figsize=(13 if wide else 9, 9.5), sharey=True, squeeze=False)
    for i, (level, col, row_title, methods) in enumerate(ROWS):
        for j, (error, error_title) in enumerate(ERRORS):
            ax = axes[i][j]
            width = 0.8 / len(methods)
            for k, m in enumerate(methods):
                name, color, marker = METHODS[m]
                xs, ys, lo, hi = [], [], [], []
                for x, s in enumerate(samples):
                    r = data[level].get((s, error, m))
                    if r is None or r[f"{col} mean"] == "":
                        continue
                    mean = float(r[f"{col} mean"])
                    sd = float(r[f"{col} sd"]) if r[f"{col} sd"] else 0.0
                    y = max(mean, FLOOR)
                    xs.append(x - 0.4 + width * (k + 0.5))
                    ys.append(y)
                    lo.append(y - max(mean - sd, FLOOR) if sd else 0)
                    hi.append(sd)
                ax.errorbar(xs, ys, yerr=[lo, hi], fmt=marker, color=color, ms=6, mec="white", mew=0.8,
                            elinewidth=1.2, capsize=0, label=name, zorder=3)
            ax.set_yscale("log")
            ax.set_ylim(FLOOR / 1.6, 3)
            ax.axhline(FLOOR, color=plot_style.MUTED, lw=0.6, ls=":", zorder=1)
            ax.set_xticks(range(len(samples)))
            ax.set_xticklabels(labels, rotation=40 if wide else 0, ha="right" if wide else "center",
                               fontsize=8 if wide else 9)
            ax.set_xlim(-0.6, len(samples) - 0.4)
            ax.grid(axis="x", visible=False)
            if i == 0:
                ax.set_title(error_title)
            if j == 0:
                ax.set_ylabel(f"{row_title}\nL1 zur Wahrheit")
            if i < len(ROWS) - 1:
                ax.tick_params(labelbottom=False)
    # Plain marker handles: errorbar handles would draw a stray whisker in the legend.
    handles = [Line2D([], [], ls="", marker=mk, color=c, ms=6, mec="white", mew=0.8) for _, c, mk in METHODS.values()]
    order = [name for name, _, _ in METHODS.values()]
    fig.legend(handles, order, loc="upper center", ncol=len(order),
               bbox_to_anchor=(0.5, 0.955), handletextpad=0.3, columnspacing=1.6)
    fig.suptitle(f"{title}: L1-Abstand zur Wahrheit (0 = perfekt, 2 = völlig daneben)", y=0.99, fontsize=11)
    fig.text(0.5, 0.005, f"Punkt = Mittelwert über die Seeds, Strich = ± eine Standardabweichung. "
             f"Logarithmische Achse; Werte unter {FLOOR:g} liegen auf der gepunkteten Linie.",
             ha="center", fontsize=8, color=plot_style.MUTED)
    fig.tight_layout(rect=(0, 0.02, 1, 0.93))
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/plots/l1_{group}.{ext}", dpi=150)
    plt.close(fig)
    print(f"written {OUT}/plots/l1_{group}.png/.pdf")
