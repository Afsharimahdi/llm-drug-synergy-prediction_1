# LLM-Drug-Synergy-Prediction

Differential perturbation encoding with a cell foundation model (Cell2Sentence)
for drug-pair synergy prediction.

## Overview

We predict whether a drug pair is synergistic or antagonistic on a given cancer
cell line. Instead of describing the cell line only by its baseline
gene-expression vector, we prompt the **C2S-Scale (Gemma-2 27B)** cell
foundation model to *generate* the cell's perturbation response to each drug and
to the pair, and we encode the difference between these generated responses
(differential perturbation encoding). A gated fusion transformer then combines
this with a Morgan-fingerprint drug encoder.

## Repository contents

| Path | What it is |
|---|---|
| `DATA (3).ipynb` | Data preparation and C2S cell-sentence generation |
| `model_4.ipynb` | Model definition, training, and the Drug+Gene vs Drug-only comparison |
| `DATA_PLOT.ipynb` | The cell-sentence rank figures |
| `c2s27b_out_0_4396_in300_out300 (2).rar` | The generated cell sentences (see *Data* below) |
| `paper_analysis/` | Scripts reproducing the descriptive numbers in the paper — see [`paper_analysis/README.md`](paper_analysis/README.md) |
| `requirements.txt` | Pinned versions used for the reported results |

## Data

- **Synergy labels:** DrugCombDB — <https://doi.org/10.1093/nar/gkz1007>
- **Gene expression:** DepMap portal, release Public 25Q2 — <https://depmap.org>
- **Generated cell sentences:** in the `.rar` archive above. Unpack it to get
  `c2s27b_out_0_4396_in300_out300 (2).csv`, which holds, for every
  (cell line, drug 1, drug 2) triplet, the untreated cell sentence and the three
  C2S-generated treated sentences.

The archive has **4,268 rows**. Applying the SMILES standardisation used for
training (largest fragment, uncharge, canonical; invalid SMILES dropped) removes
383 rows and leaves the **3,885** analysed in the paper.
`paper_analysis/recompute_filtered.py` reproduces exactly that step.

## Reproducing the paper

```bash
pip install -r requirements.txt
unrar x "c2s27b_out_0_4396_in300_out300 (2).rar"     # or any unzip tool

# descriptive results (dataset table, gene-block ranks, per-cell-line rates,
# rank-correlation decomposition)
python paper_analysis/recompute_filtered.py
python paper_analysis/spearman_all.py

# the three cell-sentence figures, at 300 dpi
python paper_analysis/regen_figs.py figures/
```

The predictive results (cross-validation and held-out test) come from
`model_4.ipynb`. Generating the cell sentences from scratch needs a GPU; we used
an RTX 5090 with the 27B model in 4-bit NF4 precision. The released archive means
you do not have to repeat that step.

## What is not in this repository yet

For transparency: the representation study reported in the paper — the one-hot
cell-line control, the identity-only and cell-line-average baselines, and the
leave-cell-out split — is not yet included here. It will be added.

## Requirements

Python 3.12; see `requirements.txt` for pinned versions. `rdkit` is needed to
reproduce the SMILES standardisation that defines the 3,885-row analysis set.

## Citation

See `CITATION.cff`. The paper is under review; the reference will be updated
once it is published.

## License

MIT — see `LICENSE`.
