import pandas as pd, numpy as np
from scipy.stats import spearmanr

SRC = r"c:/Users/afsha/Desktop/final/model_khodam/data/c2s27b_out_0_4396_in300_out300 (2).csv"
df = pd.read_csv(SRC)

i = 0
print("cell=%s  d1=%s  d2=%s  label=%s\n"
      % (df.loc[i,"cell"], df.loc[i,"drug1_name"], df.loc[i,"drug2_name"], df.loc[i,"label"]))

S = {k: df.loc[i, c].split() for k, c in
     [("Before","cell_sentence"), ("D1","c2s_treated_d1_sentence"),
      ("D2","c2s_treated_d2_sentence"), ("D12","c2s_treated_d12_direct")]}
stages = ["Before","D1","D2","D12"]

# ---- 1. how much of each sentence is even the same set of genes?
print("--- 1. gene-set overlap with Before (out of 300) ---")
for s in stages[1:]:
    print("   %-4s %3d shared, %3d new" %
          (s, len(set(S[s]) & set(S["Before"])), len(set(S[s]) - set(S["Before"]))))

# ---- 2. are the three treated sentences actually different from each other?
common = set.intersection(*[set(v) for v in S.values()])
common = sorted(common)
R = {s: np.array([S[s].index(g)+1 for g in common]) for s in stages}
print("\n--- 2. Spearman rank correlation on the %d genes present in all four ---" % len(common))
print("        " + "".join("%8s" % s for s in stages))
for a in stages:
    print("   %-5s" % a + "".join("%8.3f" % spearmanr(R[a], R[b]).statistic for b in stages))

# ---- 3. the two blocks seen in the heatmap
MT  = [g for g in common if g.startswith("MT-")]
RP  = [g for g in common if g.startswith(("RPL","RPS","RPLP"))]
print("\n--- 3. median rank of each block (lower = more expressed) ---")
print("   block            n" + "".join("%8s" % s for s in stages))
for name, gl in [("MT- (OXPHOS)", MT), ("RP  (ribosomal)", RP)]:
    idx = [common.index(g) for g in gl]
    print("   %-16s %2d" % (name, len(gl)) + "".join("%8.0f" % np.median(R[s][idx]) for s in stages))

# ---- 4. is the RP rise just displacement caused by the MT genes falling?
print("\n--- 4. controlling for displacement ---")
print("   If the %d MT genes simply drop out of the top of the list, every gene" % len(MT))
print("   below them can rise by at most that many positions.")
for s in ["D1","D2","D12"]:
    idx = [common.index(g) for g in RP]
    obs = np.median(R["Before"][idx]) - np.median(R[s][idx])
    mt_above = sum(1 for g in MT if S["Before"].index(g)+1 < np.median(R["Before"][idx]))
    print("   %-4s observed median RP rise = %3.0f positions;  MT genes above them in"
          " Before = %d  ->  displacement explains %.0f%%" % (s, obs, mt_above, 100*mt_above/max(obs,1)))

# ---- 5. how drug-specific is this, really? compare against other pairs in A2058
print("\n--- 5. same cell line (A2058), different drug pairs ---")
sub = df[df["cell"] == df.loc[i,"cell"]].head(6)
base = df.loc[i,"cell_sentence"].split()
print("   %-26s %-9s  Spearman(D12 vs Before)  median RP rank in D12" % ("drug pair","label"))
for j, row in sub.iterrows():
    d12 = row["c2s_treated_d12_direct"].split()
    com = sorted(set(d12) & set(base))
    rb = np.array([base.index(g)+1 for g in com]); rd = np.array([d12.index(g)+1 for g in com])
    rp = [k for k,g in enumerate(com) if g.startswith(("RPL","RPS","RPLP"))]
    print("   %-26s %-9s  %6.3f                  %3.0f"
          % (row["drug1_name"]+" + "+row["drug2_name"], row["label"],
             spearmanr(rb, rd).statistic, np.median(rd[rp]) if rp else -1))
