# Regenerate the three cell-sentence figures used in the paper,
# titled with the cell-line name instead of the row index.
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np, seaborn as sns, sys, os

SRC = r"c:/Users/afsha/Desktop/final/model_khodam/data/c2s27b_out_0_4396_in300_out300 (2).csv"
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)

df = pd.read_csv(SRC)
cell_id = 0
cell_name = df.loc[cell_id, "cell"]
d1_name = df.loc[cell_id, "drug1_name"]
d2_name = df.loc[cell_id, "drug2_name"]
print("cell_id=%d -> cell line %s  (%s + %s, %s)"
      % (cell_id, cell_name, d1_name, d2_name, df.loc[cell_id, "label"]))

before = df.loc[cell_id, "cell_sentence"].split()
d1  = df.loc[cell_id, "c2s_treated_d1_sentence"].split()
d2  = df.loc[cell_id, "c2s_treated_d2_sentence"].split()
d12 = df.loc[cell_id, "c2s_treated_d12_direct"].split()
print("sentence lengths:", len(before), len(d1), len(d2), len(d12))

common = set(before) & set(d1) & set(d2) & set(d12)
stages = ["Before", "D1", "D2", "D12"]
lists  = [before, d1, d2, d12]
print("genes common to all four stages:", len(common))

# ---------------------------------------------------------------- rank_all
plt.figure(figsize=(12, 8))
for g in common:
    plt.plot(stages, [l.index(g) + 1 for l in lists], marker="o", alpha=0.3, linewidth=1)
plt.xlabel("Condition", fontsize=14, fontweight="bold")
plt.ylabel("Gene Expression Rank", fontsize=14, fontweight="bold")
plt.title(f"Gene Expression Rank Changes for Cell line {cell_name}",
          fontsize=16, fontweight="bold")
plt.gca().invert_yaxis(); plt.grid(True, alpha=0.3); plt.tight_layout()
plt.savefig(f"{OUT}/rank_all.png", dpi=300); plt.close()

# ------------------------------------------------------------- rank_top20
changes = {g: [l.index(g) + 1 for l in lists] for g in common}
top = sorted(changes.items(), key=lambda x: abs(x[1][-1] - x[1][0]), reverse=True)[:20]
plt.figure(figsize=(14, 8))
for g, ranks in top:
    plt.plot(stages, ranks, marker="o", linewidth=2, label=g, alpha=0.7)
plt.xlabel("Condition", fontsize=14, fontweight="bold")
plt.ylabel("Gene Expression Rank", fontsize=14, fontweight="bold")
plt.title(f"Top 20 Genes with Largest Rank Changes (Cell line {cell_name})",
          fontsize=16, fontweight="bold")
plt.gca().invert_yaxis()
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=9)
plt.grid(True, alpha=0.3); plt.tight_layout()
plt.savefig(f"{OUT}/rank_top20.png", dpi=300); plt.close()

# ----------------------------------------------------------- rank_heatmap
rank_df = pd.DataFrame(list(changes.values()), columns=stages, index=list(changes.keys()))
rank_df["Total_Change"] = (rank_df["D12"] - rank_df["Before"]).abs()
top30 = rank_df.sort_values("Total_Change", ascending=False).head(30)
plt.figure(figsize=(10, 14))
sns.heatmap(top30[stages], cmap="RdYlGn_r", annot=True, fmt="d",
            cbar_kws={"label": "Expression Rank"}, linewidths=0.5)
plt.title(f"Top 30 Genes - Expression Rank Changes (Cell line {cell_name})",
          fontsize=14, fontweight="bold", pad=15)
plt.xlabel("Condition", fontsize=12, fontweight="bold")
plt.ylabel("Genes", fontsize=12, fontweight="bold")
plt.tight_layout()
plt.savefig(f"{OUT}/rank_heatmap.png", dpi=300); plt.close()

print("wrote rank_all.png, rank_top20.png, rank_heatmap.png ->", OUT)
