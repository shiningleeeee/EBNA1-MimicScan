#!/usr/bin/env python3
"""Figure 1: EBNA1's Gly-Ala block explains the entire apparent mimicry signal.

(a) sliding-window G+A content along EBNA1 -- residues 90-330 are 100% Gly/Ala
(b) every alignment MMseqs2 is willing to report, drawn on the query axis:
    they all sit in 350-641 and none touches the block
(c) best E-value of the real sequence vs composition-matched null families:
    the advantage over a global shuffle collapses once local order is kept

Usage: 09_plot.py  ->  reports/figure1.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

M8 = "reports/mmseqs_ebna1.m8"
PER = "reports/batch_null_per_query.csv"
REAL_BEST = 3.088


def read_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">"))


def main():
    seq = read_seq("data/cache/P03211.fasta")
    L = len(seq)
    df = pd.read_csv(M8, sep="\t", header=None,
                     names=["target", "evalue", "bits", "qstart", "qend",
                            "tstart", "tend", "alnlen", "pident", "qaln", "taln"])
    per = pd.read_csv(PER)

    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(11.5, 10.2),
        gridspec_kw={"height_ratios": [1.0, 1.35, 1.15]})

    # ---------------- (a) composition profile ----------------
    w = 20
    xs = np.arange(L - w + 1) + w / 2
    ga = np.array([sum(c in "GA" for c in seq[i:i + w]) / w for i in range(L - w + 1)])
    ax1.fill_between(xs, 0, ga, color="#2b6cb0", alpha=0.35, step="mid")
    ax1.plot(xs, ga, color="#2b6cb0", lw=0.9)
    ax1.axvspan(90, 330, color="#c53030", alpha=0.10)
    ax1.set_ylim(-0.02, 1.05)
    ax1.set_xlim(0, L)
    ax1.set_ylabel("G + A fraction\n(20-aa window)", fontsize=9)
    ax1.set_title("(a) EBNA1 is dominated by a 240-residue block that is 100% Gly/Ala "
                  "(residues 90-330)", fontsize=10.5, loc="left")
    ax1.annotate("90-330\nG+A = 1.00", xy=(210, 0.5), ha="center", fontsize=9,
                 color="#c53030", fontweight="bold")
    ax1.axhline(1.0, color="grey", lw=0.6, ls=":")
    ax1.grid(alpha=0.2)

    # ---------------- (b) where the reported alignments sit ----------------
    d = df.sort_values("evalue").reset_index(drop=True)
    ax2.axvspan(90, 330, color="#c53030", alpha=0.10)
    cmap = plt.get_cmap("viridis")
    logE = np.log10(d["evalue"].to_numpy())
    vmin, vmax = logE.min(), logE.max()
    for i, r in d.iterrows():
        c = cmap(1 - (logE[i] - vmin) / max(vmax - vmin, 1e-9))
        ax2.plot([r.qstart, r.qend], [i, i], color=c, lw=1.6, solid_capstyle="butt")
    ax2.set_xlim(0, L)
    ax2.set_ylim(-2, len(d) + 2)
    ax2.set_ylabel("alignment rank (by E-value)", fontsize=9)
    ax2.set_title(f"(b) all {len(d)} reported alignments (E {df['evalue'].min():.3g}"
                  f"-{df['evalue'].max():.2g}) sit in 350-641; none touches the block",
                  fontsize=10.5, loc="left")
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin, vmax))
    cb = fig.colorbar(sm, ax=ax2, pad=0.01, fraction=0.03)
    cb.set_label("log10 E-value (top = best)", fontsize=8)
    n_in = int(((d.qstart <= 330) & (d.qend >= 90)).sum())
    ax2.text(210, len(d) * 0.5, f"alignments overlapping 90-330: {n_in}",
             rotation=90, ha="center", va="center", fontsize=10,
             color="#c53030", fontweight="bold")

    # ---------------- (c) null families ----------------
    fams = ["w10", "w20", "w40", "w641"]
    labels = ["shuffle\nin 10-aa\nwindows", "in 20-aa\nwindows",
              "in 40-aa\nwindows", "global\nshuffle"]
    data, ps = [], []
    for f in fams:
        be = per.loc[per["family"] == f, "bestE"].to_numpy(float)
        data.append(np.clip(be, 1e-3, None))
        ps.append(float((be <= REAL_BEST).mean()))
    bp = ax3.boxplot(data, widths=0.55, showfliers=False, patch_artist=True,
                     medianprops=dict(color="black", lw=1.4))
    for patch, f in zip(bp["boxes"], fams):
        patch.set_facecolor("#90cdf4" if f != "w641" else "#fbd38d")
        patch.set_alpha(0.8)
    rng = np.random.default_rng(0)
    for i, d_ in enumerate(data):
        y = np.clip(d_, 1e-3, None)
        ax3.scatter(i + 1 + rng.normal(0, 0.045, len(y)), y, s=2.5,
                    color="#2d3748", alpha=0.25, zorder=3)
    ax3.axhline(REAL_BEST, color="#c53030", lw=2, ls="--", zorder=4)
    ax3.axhline(1.0, color="grey", lw=0.8, ls=":")
    ax3.text(4.72, 0.90, "E = 1", fontsize=8, color="grey", va="top", ha="right")
    ax3.set_yscale("log")
    ax3.set_xticks(range(1, len(fams) + 1))
    ax3.set_xticklabels([f"{l}\nP(null$\\leq$real) = {p:.3f}" for l, p in zip(labels, ps)],
                        fontsize=9)
    ax3.set_ylabel("best E-value vs 20,431\nhuman proteins", fontsize=9)
    ax3.set_title("(c) keep local composition and the advantage vanishes: "
                  "B = 1000 shuttles per family", fontsize=10.5, loc="left")
    ax3.grid(alpha=0.2, axis="y")
    ax3.set_xlim(0.4, len(fams) + 0.75)
    ax3.annotate(f"real EBNA1  E = {REAL_BEST:.2f}", xy=(len(fams) + 0.7, REAL_BEST),
                 xytext=(len(fams) + 0.7, REAL_BEST * 3.2),
                 color="#c53030", fontsize=9.5, ha="right", va="bottom",
                 fontweight="bold",
                 arrowprops=dict(arrowstyle="-", color="#c53030", lw=1.2))

    fig.suptitle("EBNA1 vs human proteome: no sequence-level mimicry beyond local composition",
                 fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig("reports/figure1.png", dpi=150)
    print("wrote reports/figure1.png")
    print(f"alignments overlapping the Gly/Ala block: {n_in}/{len(d)}")


if __name__ == "__main__":
    main()
