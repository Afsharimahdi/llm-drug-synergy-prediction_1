# Reproduce the exact filtering used in MODEL_3.ipynb / model_4.ipynb
# (4268 rows -> drop 383 invalid SMILES -> 3885) and recompute every number
# that would go into the paper, so they match Table II.
import pandas as pd, numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog("rdApp.*")

SRC = r"c:/Users/afsha/Desktop/final/model_khodam/data/c2s27b_out_0_4396_in300_out300 (2).csv"
NULLISH = {"", "none", "nan", "null", "na", "n/a", "undefined", "missing"}

def clean_smiles_string(x):
    if x is None: return None
    s = str(x).strip()
    return None if (not s) or (s.lower() in NULLISH) else s

_lfc, _unc = rdMolStandardize.LargestFragmentChooser(), rdMolStandardize.Uncharger()
def standardize(s):
    s = clean_smiles_string(s)
    if not s: return None
    mol = Chem.MolFromSmiles(s)
    if mol is None: return None
    try:
        mol = _unc.uncharge(_lfc.choose(mol))
    except Exception:
        pass
    try:
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        return None

df = pd.read_csv(SRC).dropna(subset=["label"]).reset_index(drop=True)
df["y"] = df["label"].astype(str).str.strip().str.lower().map(
    {"synergy": 1, "antagonism": 0, "1": 1, "0": 0})
print("Shape raw:", df.shape, "| cells:", df["cell"].nunique(),
      "| label:", df["y"].value_counts().to_dict())

s1 = df["drug1_smiles"].map(standardize)
s2 = df["drug2_smiles"].map(standardize)
ok = s1.notna() & s2.notna()
print("SMILES standardization: dropped %d invalid rows" % (~ok).sum())
df = df[ok].copy()
df["drug1_smiles"], df["drug2_smiles"] = s1[ok], s2[ok]
df = df.reset_index(drop=True)
print("After standardize:", df.shape, "-> matches paper's 3,885:", len(df) == 3885)
print("synergy rate: %.3f   (paper Table II: 0.319)" % df["y"].mean())
print("cell lines:   %d      (paper Table II: 52)" % df["cell"].nunique())
n_drugs = len(set(df["drug1_smiles"]) | set(df["drug2_smiles"]))
print("distinct drugs: %d    (paper Discussion: 106)" % n_drugs)

# ============================================================ item 3
print("\n" + "=" * 62)
print("ITEM 3 -- synergy rate per cell line (on the filtered 3,885)")
g = df.groupby("cell")["y"].agg(n="size", rate="mean")
g20 = g[g.n >= 20].sort_values("rate")
print("  overall rate            : %.3f" % df["y"].mean())
print("  cell lines with n>=20   : %d  (of %d)" % (len(g20), len(g)))
print("  min / median / max rate : %.3f / %.3f / %.3f"
      % (g20.rate.min(), g20.rate.median(), g20.rate.max()))
print("\n  lowest 3:\n" + g20.head(3).round(3).to_string())
print("\n  highest 3:\n" + g20.tail(3).round(3).to_string())
for c in ["A2058", "A2780", "A427"]:
    if c in g.index:
        print("  %-7s n=%4d  rate=%.3f  (%.1f%% of the dataset)"
              % (c, g.loc[c, "n"], g.loc[c, "rate"], 100 * g.loc[c, "n"] / len(df)))

# ============================================================ items 1 & 2
print("\n" + "=" * 62)
row = df[(df["cell"] == "A2058") & (df["drug1_name"] == "5-FU")
         & (df["drug2_name"] == "ABT-888")]
if len(row) == 0:
    print("the illustrated triplet did NOT survive SMILES filtering"); raise SystemExit
row = row.iloc[0]
print("ITEMS 1 & 2 -- %s : %s + %s (%s), survives filtering"
      % (row["cell"], row["drug1_name"], row["drug2_name"], row["label"]))

S = {k: row[c].split() for k, c in
     [("Before", "cell_sentence"), ("D1", "c2s_treated_d1_sentence"),
      ("D2", "c2s_treated_d2_sentence"), ("D12", "c2s_treated_d12_direct")]}
stages = ["Before", "D1", "D2", "D12"]
common = sorted(set.intersection(*[set(v) for v in S.values()]))
R = {s: np.array([S[s].index(g) + 1 for g in common]) for s in stages}

print("\n  median rank by block (n common genes = %d)" % len(common))
print("  %-26s %3s" % ("block", "n") + "".join("%8s" % s for s in stages))
for name, pref in [("Mitochondrial (MT-)", ("MT-",)),
                   ("Ribosomal protein", ("RPL", "RPS", "RPLP"))]:
    idx = [i for i, g in enumerate(common) if g.startswith(pref)]
    print("  %-26s %3d" % (name, len(idx))
          + "".join("%8.0f" % np.median(R[s][idx]) for s in stages))

print("\n  Spearman rank correlation")
print("  %-8s" % "" + "".join("%8s" % s for s in stages))
for a in stages:
    print("  %-8s" % a + "".join("%8.3f" % spearmanr(R[a], R[b]).statistic for b in stages))

# ============================================================ the dead idea
print("\n" + "=" * 62)
print("REJECTED IDEA -- deviation from baseline as a predictor (filtered set)")
dev = []
for _, r in df.iterrows():
    b, d = r["cell_sentence"].split(), r["c2s_treated_d12_direct"].split()
    pd_ = {g: i for i, g in enumerate(d)}
    com = [g for g in b if g in pd_]
    if len(com) < 20: dev.append(np.nan); continue
    rb = np.array([b.index(g) for g in com]); rd = np.array([pd_[g] for g in com])
    dev.append(1 - spearmanr(rb, rd).statistic)
df["dev"] = dev
m = df["dev"].notna()
print("  ROC-AUC of deviation alone = %.4f   (0.50 = chance)"
      % roc_auc_score(df.loc[m, "y"], df.loc[m, "dev"]))
print("  -> still uninformative; keep it out of the paper.")
