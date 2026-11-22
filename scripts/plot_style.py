"""Shared look of the evaluation figures: one fixed colour and marker per method, in the order the methods
are listed (naive renorm., naive, then the topic models as given), so a method keeps its colour in every
figure. Colours are the first slots of a categorical palette checked for colour-vision deficiency; since two
of them are light against white, every method also has its own marker shape and a legend, so colour never
carries identity alone."""
import matplotlib

matplotlib.use("Agg")  # files only, no window
import matplotlib.pyplot as plt  # noqa: E402

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MARKERS = ["o", "s", "^", "D", "v", "p", "h", "<"]  # solid shapes: thin ones (P, X, *) vanish under the white edge
INK = "#1f1f1e"  # text
MUTED = "#6b6a63"  # axes, ticks, reference lines
GRID = "#e4e3dd"


def style(method_index):
    i = method_index % len(COLORS)
    return {"color": COLORS[i], "marker": MARKERS[i]}


def setup():
    plt.rcParams.update({
        "font.size": 9,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "axes.titlesize": 10,
        "axes.titlecolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "pdf.fonttype": 42,  # embeds text as text, so it stays selectable in the PDF
    })
    return plt
