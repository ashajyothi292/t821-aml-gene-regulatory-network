"""
gsea_preranked.py

Pre-ranked gene set enrichment analysis of the three RNA-seq comparisons
(all patients, core, other vs control) against Reactome pathways, followed by a
dot plot of the significant pathways (Figure 14 of the report).

Input:  rnaseq/deseq2/ranks_<comparison>.rnk  (from deseq2_rnaseq.R)
        a Reactome GMT file from MSigDB (c2.cp.reactome.*.symbols.gmt)
Output: rnaseq/gsea/<comparison>/  full GSEApy results
        rnaseq/gsea/reactome_significant.csv
        rnaseq/gsea/reactome_dotplot.png

Usage:  python scripts/06_rnaseq/gsea_preranked.py <project_dir> <reactome.gmt>
"""

import sys
from pathlib import Path

import gseapy as gp
import matplotlib.pyplot as plt
import pandas as pd

FDR_CUTOFF = 0.25
COMPARISONS = {"all": "Control vs t821 all",
               "core": "Control vs t821 core",
               "other": "Control vs t821 other"}


def main(project_dir: Path, gmt: Path) -> None:
    deseq_dir = project_dir / "rnaseq" / "deseq2"
    out_dir = project_dir / "rnaseq" / "gsea"
    out_dir.mkdir(parents=True, exist_ok=True)

    tables = []
    for key, label in COMPARISONS.items():
        ranks = pd.read_csv(deseq_dir / f"ranks_{key}.rnk", sep="\t", header=None,
                            names=["gene", "score"])
        res = gp.prerank(rnk=ranks, gene_sets=str(gmt), permutation_num=1000,
                         min_size=15, max_size=500, seed=42, threads=4,
                         outdir=str(out_dir / key), verbose=False)
        df = res.res2d.copy()
        df["comparison"] = label
        tables.append(df)
        n_sig = (df["FDR q-val"].astype(float) < FDR_CUTOFF).sum()
        print(f"{label}: {len(df)} pathways tested, {n_sig} with FDR < {FDR_CUTOFF}")

    allres = pd.concat(tables, ignore_index=True)
    allres["FDR"] = allres["FDR q-val"].astype(float)
    allres["NES"] = allres["NES"].astype(float)
    sig = allres[allres["FDR"] < FDR_CUTOFF].copy()
    sig["Pathway"] = (sig["Term"].str.replace("REACTOME_", "", regex=False)
                      .str.replace("_", " ").str.capitalize())
    sig.to_csv(out_dir / "reactome_significant.csv", index=False)

    if sig.empty:
        print("No pathway passed the FDR cut-off; no plot drawn.")
        return

    # dot plot: one column per comparison, red = up in AML, blue = down, size = |NES|
    pathways = sig.groupby("Pathway")["NES"].mean().sort_values().index.tolist()
    xpos = {label: i for i, label in enumerate(COMPARISONS.values())}
    fig, ax = plt.subplots(figsize=(9, max(4, 0.35 * len(pathways))))
    ax.scatter(sig["comparison"].map(xpos), sig["Pathway"].map(pathways.index),
               s=(sig["NES"].abs() ** 2) * 60,
               c=["#E74C3C" if v > 0 else "#2980B9" for v in sig["NES"]],
               edgecolors="white", linewidths=0.6)
    ax.set_xticks(list(xpos.values()), list(xpos.keys()))
    ax.set_yticks(range(len(pathways)), pathways, fontsize=8)
    ax.set_xlim(-0.6, len(xpos) - 0.4)
    ax.set_title(f"Enriched Reactome pathways (FDR < {FDR_CUTOFF})")
    ax.grid(axis="y", color="#eeeeee")
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out_dir / "reactome_dotplot.png", dpi=300)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(Path(sys.argv[1]), Path(sys.argv[2]))
