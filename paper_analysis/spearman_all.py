# Does the "shared component dominates" pattern hold across the whole dataset,
# or only in the illustrated triplet?  Run on the paper's filtered 3,885 rows.
import pandas as pd, numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from scipy.stats import spearmanr, wilcoxon
RDLogger.DisableLog("rdApp.*")

SRC = r"c:/Users/afsha/Desktop/final/model_khodam/data/c2s27b_out_0_4396_in300_out300 (2).csv"
NULLISH = {"", "none", "nan", "null", "na", "n/a", "undefined", "missing"}
_lfc, _unc = rdMolStandardize.LargestFragmentChooser(), rdMolStandardize.Uncharger()

def standardize(x):
    if x is None: return None
    s = str(x).strip()
    if (not s) or (s.lower() in NULLISH): return None
    mol = Chem.MolFromSmiles(s)
    if mol is None: return None
    try: mol = _unc.uncharge(_lfc.choose(mol))
    except Exception: pass
    try: return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception: return None

df = pd.read_csv(SRC).dropna(subset=["label"]).reset_index(drop=True)
ok = df["drug1_smiles"].map(standardize).notna() & df["drug2_smiles"].map(standardize).notna()
df = df[ok].reset_index(drop=True)
print("filtered rows:", len(df))

cols = ["cell_sentence", "c2s_treated_d1_sentence",
        "c2s_treated_d2_sentence", "c2s_treated_d12_direct"]
base_treated, treated_treated, n_common = [], [], []

for _, r in df.iterrows():
    S = [str(r[c]).split() for c in cols]
    if any(len(s) < 50 for s in S):
        base_treated.append(np.nan); treated_treated.append(np.nan); n_common.append(0); continue
    common = sorted(set(S[0]) & set(S[1]) & set(S[2]) & set(S[3]))
    if len(common) < 30:
        base_treated.append(np.nan); treated_treated.append(np.nan); n_common.append(len(common)); continue
    pos = [{g: i for i, g in enumerate(s)} for s in S]
    R = [np.array([p[g] for g in common]) for p in pos]
    bt = [spearmanr(R[0], R[k]).statistic for k in (1, 2, 3)]           # baseline vs treated
    tt = [spearmanr(R[a], R[b]).statistic for a, b in ((1,2),(1,3),(2,3))]  # treated vs treated
    base_treated.append(np.mean(bt)); treated_treated.append(np.mean(tt)); n_common.append(len(common))

df["bt"] = base_treated
df["tt"] = treated_treated
df["n_common"] = n_common
d = df.dropna(subset=["bt", "tt"])
print("rows with a usable comparison:", len(d))
print("genes common to all four sentences: median %d  (IQR %d-%d)"
      % (d.n_common.median(), d.n_common.quantile(.25), d.n_common.quantile(.75)))

print("\n=== mean Spearman across the dataset ===")
print("  baseline vs treated   : %.3f  (sd %.3f, IQR %.3f-%.3f)"
      % (d.bt.mean(), d.bt.std(), d.bt.quantile(.25), d.bt.quantile(.75)))
print("  treated  vs treated   : %.3f  (sd %.3f, IQR %.3f-%.3f)"
      % (d.tt.mean(), d.tt.std(), d.tt.quantile(.25), d.tt.quantile(.75)))
print("  difference (tt - bt)  : %.3f" % (d.tt - d.bt).mean())

frac = (d.tt > d.bt).mean()
print("\n  rows where treated-treated > baseline-treated: %.1f%%  (%d of %d)"
      % (100 * frac, (d.tt > d.bt).sum(), len(d)))
w = wilcoxon(d.tt, d.bt)
print("  Wilcoxon signed-rank: statistic=%.4g, p=%.3g" % (w.statistic, w.pvalue))

print("\n=== the illustrated triplet, for comparison ===")
t = df[(df.cell == "A2058") & (df.drug1_name == "5-FU") & (df.drug2_name == "ABT-888")]
if len(t):
    t = t.iloc[0]
    print("  A2058 / 5-FU + ABT-888: baseline-treated %.3f, treated-treated %.3f"
          % (t.bt, t.tt))
    pct = (d.bt < t.bt).mean() * 100
    print("  its baseline-treated value sits at the %.0fth percentile of the dataset" % pct)
