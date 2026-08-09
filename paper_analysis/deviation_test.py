# Hypothesis: "how far the generated combination response departs from the
# untreated baseline" is by itself predictive of the synergy label.
import pandas as pd, numpy as np
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.metrics import roc_auc_score

SRC = r"c:/Users/afsha/Desktop/final/model_khodam/data/c2s27b_out_0_4396_in300_out300 (2).csv"
df = pd.read_csv(SRC).dropna(
    subset=["cell_sentence", "c2s_treated_d12_direct", "label"]).reset_index(drop=True)
print("rows:", len(df), "| labels:", dict(df["label"].value_counts()))

rho, n_new, med_shift = [], [], []
for _, r in df.iterrows():
    b = r["cell_sentence"].split()
    d = r["c2s_treated_d12_direct"].split()
    pos_b = {g: i for i, g in enumerate(b)}
    pos_d = {g: i for i, g in enumerate(d)}
    common = [g for g in b if g in pos_d]
    if len(common) < 20:
        rho.append(np.nan); n_new.append(np.nan); med_shift.append(np.nan); continue
    rb = np.array([pos_b[g] for g in common])
    rd = np.array([pos_d[g] for g in common])
    rho.append(spearmanr(rb, rd).statistic)
    n_new.append(len(set(d) - set(b)))
    med_shift.append(np.median(np.abs(rd - rb)))

df["rho"] = rho                       # 1.0 = response identical to baseline
df["deviation"] = 1 - df["rho"]       # bigger = the drugs changed the cell more
df["n_new_genes"] = n_new
df["median_abs_shift"] = med_shift
df = df.dropna(subset=["rho"])
y = (df["label"] == "synergy").astype(int).values

print("\n--- deviation from baseline, by label ---")
print(df.groupby("label")[["deviation", "n_new_genes", "median_abs_shift"]]
        .agg(["mean", "median", "std", "count"]).round(4).to_string())

print("\n--- is each single feature predictive on its own? ---")
print("   %-18s %8s %8s %10s" % ("feature", "AUC", "U-p", "direction"))
for f in ["deviation", "n_new_genes", "median_abs_shift"]:
    v = df[f].values
    auc = roc_auc_score(y, v)
    p = mannwhitneyu(v[y == 1], v[y == 0]).pvalue
    print("   %-18s %8.4f %8.2e %10s"
          % (f, auc, p, "higher->synergy" if auc > .5 else "lower->synergy"))

print("\n(AUC 0.50 = no information. The paper's full model reaches ~0.73-0.77.)")

# the same check restricted to A2058, the cell line in the figure
a = df[df["cell"] == "A2058"]
if len(a) > 20 and a["label"].nunique() > 1:
    ya = (a["label"] == "synergy").astype(int).values
    print("\n--- A2058 only (n=%d, %d synergy / %d antagonism) ---"
          % (len(a), ya.sum(), len(ya) - ya.sum()))
    print("   deviation AUC = %.4f" % roc_auc_score(ya, a["deviation"].values))
    print(a.groupby("label")["deviation"].agg(["mean", "median", "count"]).round(4).to_string())
