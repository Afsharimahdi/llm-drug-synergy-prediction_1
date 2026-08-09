# LLM-Drug-Synergy-Prediction

Differential perturbation encoding with a cell foundation model (Cell2Sentence)
for drug-pair synergy prediction.

## Overview
We predict whether a drug pair is synergistic or antagonistic on a given cancer
cell line. Instead of describing the cell line only by its baseline gene-expression
vector, we prompt the **C2S-Scale (Gemma-2 27B)** cell foundation model to *generate*
the cell's perturbation response to each drug, and we encode the difference between
these generated responses (differential perturbation encoding). A gated fusion
transformer then combines this with a Morgan-fingerprint drug encoder.

## Repository contents
- `DATA (3).ipynb` — data preparation and C2S cell-sentence generation
- `DATA_PLOT.ipynb` — the cell-sentence rank figures
- `model_4.ipynb` — model definition, training and evaluation
- `c2s27b_out_0_4396_in300_out300 (2).rar` — the generated cell sentences
- `paper_analysis/` — scripts reproducing the descriptive numbers in the paper
  (see `paper_analysis/README.md`)

## Data
- Synergy labels: DrugCombDB (https://doi.org/10.1093/nar/gkz1007)
- Gene expression: DepMap portal, release Public 25Q2 (https://depmap.org)
- C2S-generated cell sentences are produced by the data notebook.

## Requirements
Python 3.12, PyTorch 2.10, RDKit, scikit-learn, NumPy, pandas, transformers,
bitsandbytes, accelerate. (Generating C2S sentences needs a GPU; we used an RTX 5090.)

## Citation
If you use this code, please cite the paper (see the repository or contact the author).

## License
MIT
