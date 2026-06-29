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
- `DATA__3_.ipynb` — data preparation and C2S cell-sentence generation
- `model_4.ipynb` — model definition and training
- `leave_cell_out_eval.py` — leave-cell-out evaluation (Drug+Gene vs Drug-only)
- `multiseed_delong_eval.py` — multi-seed + DeLong significance evaluation

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
