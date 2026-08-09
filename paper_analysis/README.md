# paper_analysis

Scripts that reproduce the descriptive numbers reported in the paper. They are
kept separate from the model code because none of them train anything: they read
the generated cell sentences and recompute statistics.

All of them expect the C2S output table (unpack
`c2s27b_out_0_4396_in300_out300 (2).rar` in the repository root) and reproduce
the same SMILES standardisation used for training, so the row counts match the
dataset table in the paper: **4,268 raw rows → 383 dropped for invalid SMILES →
3,885 analysed**.

## What each script produces

| Script | Reproduces |
|---|---|
| `recompute_filtered.py` | The dataset filter (4,268 → 3,885), the median gene-block ranks reported in the gene-block table, the Spearman correlations for the illustrated triplet, and the per-cell-line synergy rates quoted in the representation-study section |
| `spearman_all.py` | The dataset-wide rank-correlation decomposition quoted in the discussion: mean Spearman 0.908 (baseline vs treated) against 0.945 (treated vs treated), holding in 84.6% of triplets |
| `regen_figs.py` | The three cell-sentence figures (`rank_all`, `rank_top20`, `rank_heatmap`) at 300 dpi |
| `bio_check.py` | The exploratory block analysis that preceded the gene-block table |
| `deviation_test.py` | A **negative** result that is deliberately not in the paper: using the magnitude of the generated response as a single predictor of synergy gives ROC-AUC 0.475, i.e. no signal |
| `verify_refs.py` | Cross-checks the bibliography against local PDFs of the cited papers (run it from the folder holding those PDFs) |

## Running them

```bash
python recompute_filtered.py
python spearman_all.py
python regen_figs.py <output-directory>
```

Each script has the path to the C2S output table at the top; edit it to match
where you unpacked the archive.

## Requirements

```
python >= 3.12
pandas  scipy  numpy  scikit-learn  rdkit  matplotlib  seaborn
```

`rdkit` is needed only to reproduce the SMILES standardisation that defines the
3,885-row analysis set.
