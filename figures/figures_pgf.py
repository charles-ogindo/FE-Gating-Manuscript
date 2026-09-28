"""Publication figures for the convergence-gating manuscript, colour version.

Figure 1  Ranking performance across the fourteen benchmark receptors.
Figure 2  Effect of the retention threshold on the reported correlation.
Figure 3  Table of contents graphic.

Run this file. It writes six files into the working directory, a vector PDF
and a 600 dpi PNG for each figure.
"""
import matplotlib
matplotlib.use("pgf")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#33373d",
    "text.color": "#1d2125",
    "axes.labelcolor": "#1d2125",
    "xtick.color": "#33373d",
    "ytick.color": "#33373d",
    "savefig.dpi": 600,
    "pgf.texsystem": "pdflatex",
    "pgf.rcfonts": False,
    "pgf.preamble": r"\usepackage[T1]{fontenc}",
    "savefig.bbox": "tight",
})

TEAL = "#0d7d87"
BLUE = "#1f5fa9"
GOLD = "#e0a200"
CORAL = "#d1495b"
GREEN = "#2e8b57"
SLATE = "#4a5568"
PAPER = "#f4f1ea"

# ----------------------------------------------------------------------
# Figure 1.  Ranking performance, fourteen receptors
# ----------------------------------------------------------------------
# receptor, n, rho gated, p gated, declined
DATA = [
    ("cMET",      5, +0.900, 0.037,  0),
    ("TYK2",     13, +0.780, 0.002,  0),
    ("MCL1",     25, +0.745, 0.00005, 2),
    ("PTP1B",    22, +0.684, 0.001,  3),
    ("Thrombin", 22, +0.619, 0.003,  1),
    ("PFKFB3",   32, +0.561, 0.0008, 0),
    ("CDK8",     31, +0.545, 0.002,  0),
    ("SHP2",     24, +0.521, 0.009,  0),
    ("HIF2A",    37, +0.461, 0.004,  0),
    ("EG5",      27, +0.371, 0.057,  0),
    ("TNKS2",    27, +0.366, 0.066,  1),
    ("p38",      29, +0.330, 0.081,  0),
    ("SYK",      44, +0.219, 0.163,  2),
    ("PDE2",     21, -0.341, 0.141,  1),
]

cmap = LinearSegmentedColormap.from_list(
    "rho", [CORAL, "#e8b4ac", "#dfe5e8", "#7fb8c4", TEAL, "#075a63"])
norm = Normalize(vmin=-0.45, vmax=0.95)

fig, ax = plt.subplots(figsize=(5.0, 4.4))
ax.set_facecolor("#fbfaf7")

names = [d[0] for d in DATA]
rho = np.array([d[2] for d in DATA])
pv = np.array([d[3] for d in DATA])
nn = [d[1] for d in DATA]
dec = [d[4] for d in DATA]
y = np.arange(len(DATA))[::-1]

ax.axvspan(-0.55, 0.0, color="#f6e7e9", zorder=0)
ax.axvline(0.0, color=SLATE, lw=1.0, zorder=2)

for yi, r in zip(y, rho):
    ax.plot([0, r], [yi, yi], color=cmap(norm(r)), lw=2.6, alpha=0.55,
            solid_capstyle="round", zorder=3)

sig = pv < 0.05
ax.scatter(rho[sig], y[sig], s=78, c=[cmap(norm(v)) for v in rho[sig]],
           edgecolors="white", linewidths=1.3, zorder=5)
ax.scatter(rho[~sig], y[~sig], s=78, facecolors="white",
           edgecolors=[cmap(norm(v)) for v in rho[~sig]], linewidths=1.8,
           zorder=5)

for yi, r, d in zip(y, rho, dec):
    if d:
        ax.annotate("%d declined" % d, (r, yi), textcoords="offset points",
                    xytext=(12, -0.5), ha="left", va="center",
                    fontsize=6.5, color=SLATE, style="italic")

ax.set_yticks(y)
ax.set_yticklabels(["%s  (%d)" % (nm, n) for nm, n in zip(names, nn)],
                   fontsize=8.5)
ax.set_xlabel("Spearman correlation with measured affinity", fontsize=9)
ax.set_xlim(-0.55, 1.10)
ax.set_ylim(-0.9, len(DATA) - 0.1)
for side in ("top", "right", "left"):
    ax.spines[side].set_visible(False)
ax.tick_params(length=3)
ax.grid(axis="x", color="#dfe3e6", lw=0.6, zorder=1)
ax.set_axisbelow(True)

filled = plt.Line2D([], [], marker="o", ls="", ms=7, mfc=TEAL,
                    mec="white", mew=1.2, label="$P < 0.05$")
hollow = plt.Line2D([], [], marker="o", ls="", ms=7, mfc="white",
                    mec=SLATE, mew=1.5, label="$P \\geq 0.05$")
ax.legend(handles=[filled, hollow], frameon=False, loc="lower right",
          fontsize=8, handletextpad=0.5)

fig.savefig("figure1_ranking.pgf")
plt.close(fig)

# ----------------------------------------------------------------------
# Figure 2.  Retention threshold
# ----------------------------------------------------------------------
SWEEP = {
    "cMET":     {0.60: (+0.900, 0),  0.70: (+0.900, 0),  0.80: (None,   1)},
    "TYK2":     {0.60: (+0.780, 0),  0.70: (+0.780, 0),  0.80: (+0.780, 0)},
    "MCL1":     {0.60: (+0.738, 0),  0.70: (+0.745, 2),  0.80: (+0.462, 12)},
    "PTP1B":    {0.60: (+0.684, 3),  0.70: (+0.684, 3),  0.80: (+0.684, 3)},
    "Thrombin": {0.60: (+0.638, 0),  0.70: (+0.619, 1),  0.80: (+0.548, 14)},
    "PFKFB3":   {0.60: (+0.561, 0),  0.70: (+0.561, 0),  0.80: (+0.558, 1)},
    "CDK8":     {0.60: (+0.545, 0),  0.70: (+0.545, 0),  0.80: (+0.508, 1)},
    "SHP2":     {0.60: (+0.521, 0),  0.70: (+0.521, 0),  0.80: (+0.521, 0)},
    "HIF2A":    {0.60: (+0.461, 0),  0.70: (+0.461, 0),  0.80: (+0.461, 0)},
    "TNKS2":    {0.60: (+0.366, 1),  0.70: (+0.366, 1),  0.80: (+0.230, 3)},
    "EG5":      {0.60: (+0.371, 0),  0.70: (+0.371, 0),  0.80: (+0.371, 0)},
    "p38":      {0.60: (+0.330, 0),  0.70: (+0.330, 0),  0.80: (+0.330, 0)},
    "SYK":      {0.60: (+0.219, 2),  0.70: (+0.219, 2),  0.80: (+0.219, 2)},
    "PDE2":     {0.60: (-0.336, 0),  0.70: (-0.341, 1),  0.80: (-0.364, 6)},
}
TH = [0.60, 0.70, 0.80]
MOVERS = {"MCL1": CORAL, "Thrombin": "#7b5ea7", "TNKS2": GOLD,
          "CDK8": BLUE, "PDE2": "#b5476a", "PFKFB3": TEAL}
OFF1 = {"PFKFB3": 10, "Thrombin": 0, "CDK8": -10, "MCL1": -20, "TNKS2": -2, "PDE2": -2}
OFF2 = {"CDK8": 6, "PFKFB3": -8, "Thrombin": -2, "MCL1": -2, "TNKS2": -2, "PDE2": -2}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.3, 3.1))

for a in (ax1, ax2):
    a.set_facecolor("#fbfaf7")
    a.axvspan(0.585, 0.715, color="#e7f0f1", zorder=0)
    a.set_xticks(TH)
    a.set_xlim(0.565, 0.885)
    a.set_xlabel("Retention threshold", fontsize=9)
    for side in ("top", "right"):
        a.spines[side].set_visible(False)
    a.grid(axis="y", color="#e4e8ea", lw=0.6, zorder=1)
    a.set_axisbelow(True)
    a.tick_params(length=3)

for a in (ax1, ax2):
    a.set_facecolor("#fbfaf7")
    a.axvspan(0.585, 0.715, color="#e7f0f1", zorder=0)
    a.set_xticks(TH)
    a.set_xlim(0.568, 0.895)
    a.set_xlabel("Retention threshold", fontsize=9)
    for side in ("top", "right"):
        a.spines[side].set_visible(False)
    a.grid(axis="y", color="#e4e8ea", lw=0.6, zorder=1)
    a.set_axisbelow(True)
    a.tick_params(length=3)

for name, rows in SWEEP.items():
    xs = [t for t in TH if rows[t][0] is not None]
    ys = [rows[t][0] for t in xs]
    if name in MOVERS:
        c = MOVERS[name]
        ax1.plot(xs, ys, marker="o", ms=5, lw=1.9, color=c, zorder=5,
                 markeredgecolor="white", markeredgewidth=0.9)
        ax1.annotate(name, (xs[-1], ys[-1]), textcoords="offset points",
                     xytext=(6, OFF1.get(name, -2)), fontsize=7.5, color=c,
                     fontweight="bold")
    else:
        ax1.plot(xs, ys, marker="o", ms=3.4, lw=1.1, color="#b9c2c8", zorder=3)
ax1.axhline(0.0, color=SLATE, lw=0.8, ls=":", zorder=2)
ax1.set_ylabel("Spearman correlation", fontsize=9)
ax1.set_title("Ranking is unchanged from 0.60 to 0.70", fontsize=8.5,
              color=SLATE, pad=7)

for name, rows in SWEEP.items():
    ys = [rows[t][1] for t in TH]
    if name in MOVERS:
        c = MOVERS[name]
        ax2.plot(TH, ys, marker="o", ms=5, lw=1.9, color=c, zorder=5,
                 markeredgecolor="white", markeredgewidth=0.9)
        ax2.annotate(name, (TH[-1], ys[-1]), textcoords="offset points",
                     xytext=(6, OFF2.get(name, -2)), fontsize=7.5, color=c,
                     fontweight="bold")
    else:
        ax2.plot(TH, ys, marker="o", ms=3.4, lw=1.1, color="#b9c2c8", zorder=3)
ax2.set_ylabel("Trajectories declined", fontsize=9)
ax2.set_ylim(-1.2, 15.6)
ax2.set_title("0.80 removes trajectories that carry signal", fontsize=8.5,
              color=SLATE, pad=7)
grey = plt.Line2D([], [], color="#b9c2c8", lw=1.1, marker="o", ms=3.4,
                  label="unchanged at 0.80")
ax2.legend(handles=[grey], frameon=False, loc="upper left", fontsize=7.5,
           handletextpad=0.6)

fig.tight_layout()
fig.savefig("figure2_threshold.pgf")
plt.close(fig)

# ----------------------------------------------------------------------
# Figure 3.  Table of contents graphic
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(3.6, 1.9))
ax.set_xlim(-0.2, 10.8)
ax.set_ylim(0, 5)
ax.axis("off")
fig.patch.set_facecolor("white")


def box(x, y, w, h, text, fc, ec, fs=8, weight="normal", tc="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.10,rounding_size=0.14",
                                linewidth=1.2, edgecolor=ec, facecolor=fc,
                                zorder=3))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, fontweight=weight, color=tc, zorder=4)


def arrow(x1, y1, x2, y2, c=SLATE, lw=1.6):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=11, linewidth=lw, color=c,
                                 zorder=2))


box(0.10, 2.10, 1.95, 0.95, "Docking", "#dce8f4", BLUE, fs=8, tc=BLUE,
    weight="bold")
arrow(2.25, 2.57, 2.78, 2.57, BLUE)
box(2.95, 2.10, 1.70, 0.95, "MD", "#dce8f4", BLUE, fs=8, tc=BLUE,
    weight="bold")
arrow(4.85, 2.57, 5.35, 2.57, BLUE)

ax.add_patch(FancyBboxPatch((5.40, 1.32), 2.50, 2.50,
                            boxstyle="round,pad=0.10,rounding_size=0.16",
                            linewidth=2.0, edgecolor=GOLD,
                            facecolor="#fdf3d8", zorder=3))
ax.text(6.65, 3.38, "Gate", ha="center", va="center", fontsize=9.5,
        fontweight="bold", color="#9a6f00", zorder=4)
ax.text(6.65, 2.62, "pose RMSD\nconverged", ha="center", va="center",
        fontsize=6.6, color="#6b4e00", zorder=4)
ax.text(6.65, 1.80, "core contacts\nheld", ha="center", va="center",
        fontsize=6.6, color="#6b4e00", zorder=4)

arrow(8.02, 3.20, 8.52, 3.50, GREEN)
box(8.62, 3.15, 1.55, 0.78, "$\\Delta G$", "#dff0e5", GREEN, fs=9,
    tc=GREEN, weight="bold")
arrow(8.02, 1.95, 8.52, 1.62, CORAL)
box(8.62, 1.05, 1.55, 0.78, "declined", "#fbe4e7", CORAL, fs=6.5, tc=CORAL,
    weight="bold")

ax.text(5.2, 0.40, "359 compounds, 14 receptors", ha="center", va="center",
        fontsize=7.5, style="italic", color=SLATE)

fig.savefig("figure3_toc.pgf")
plt.close(fig)

print("figures written")
