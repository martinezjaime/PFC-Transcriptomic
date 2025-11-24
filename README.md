# PFC-Transcriptomic Clocks Repository

## Overview

This repository provides a comprehensive and reproducible framework for constructing, evaluating, and benchmarking transcriptomic aging clocks using gene expression data from the prefrontal cortex (PFC). It integrates multiple aging clock methodologies, including RNAAgeCalc, stochastic epigenetic clock simulations, epigenetic-transcriptomic correlation analyses, and retraining of transcriptomic clocks under different modeling strategies.

## Key Features

* **RNAAgeCalc Implementation**: End-to-end workflow for running RNAAgeCalc on bulk RNA-seq datasets, including preprocessing, normalization, age prediction, and performance evaluation.
* **Stochastic Clock Simulations**: Custom simulation framework implementing stochastic aging mechanisms (StocH, StocZ, StocP) adapted from epigenetic clock models to transcriptomic aging.
* **Epigenetic–Transcriptomic Correlation Analysis**: Tools to quantify and visualize correlations between DNA methylation age markers and gene expression signatures across PFC datasets.
* **Retraining of Transcriptomic Clocks**: Modular pipelines to retrain or fine-tune transcriptomic aging models using elastic net regression, deep learning, or pathway-informed architectures.
* **Reproducible, Modular Code Structure**.

## Directory Structure

```
├── data/                 # Input RNA-seq or processed expression matrices
├── notebooks/            # Jupyter/RMarkdown notebooks for exploratory analysis
├── scripts/              # Core analysis scripts (Python/R)
├── models/               # Saved models, coefficients, and retrained clocks
├── simulations/          # Stochastic clock simulation utilities
├── results/              # Output tables, metrics, and figures
└── README.md             # Project documentation
```

## Workflows

### 1. RNAAgeCalc Pipeline

* Data preprocessing (QC, filtering, normalization)
* Age prediction with RNAAgeCalc
* Model performance assessment (MAE, RMSE, correlation)
* Residualization and downstream association testing

### 2. Stochastic Clock Simulations

* Simulation of age-related trajectories based on defined stochastic rules
* Generation of synthetic expression matrices
* Elastic net training (α = 0.5) and cross-validation
* Stability analysis of clock coefficients

### 3. Epigenetic–Transcriptomic Correlations

* Integration of CpG-based epigenetic age scores
* Correlation and partial correlation analyses
* Gene set and pathway-level enrichments

### 4. Retraining and Model Extension

* Baseline linear models (EN, LASSO, Ridge)
* Deep learning architectures (TF/Keras, pathway-informed models)
* Evaluation on independent datasets
* Interpretability modules

## Installation

```bash
git clone https://github.com/yourusername/pfc-transcriptomic-clocks.git
cd pfc-transcriptomic-clocks
```

### Dependencies

* Python ≥ 3.10
* R ≥ 4.3
* Required Python packages: `pandas`, `numpy`, `scikit-learn`, `tensorflow`, `matplotlib`, `seaborn`
* Required R packages: `RNAAgeCalc`, `glmnet`, `Seurat`, `tidyverse`

## Usage

### Running RNAAgeCalc

```
Rscript scripts/run_rnaagecalc.R --input data/expression.tsv --metadata data/pheno.tsv
```

### Running Stochastic Simulations

```
python scripts/run_stochastic_clock.py --n_sim 1000 --model StocH
```

### Retraining a Transcriptomic Clock

```
python scripts/retrain_clock.py --input data/expression.tsv --age data/age.tsv
```

## Outputs

* **Age predictions** for each method
* **Simulation-derived coefficients**
* **Correlation matrices**
* **Figures** (MAE comparisons, coefficient stability, correlation heatmaps)

## Citation

If you use this repository, please cite:

```
Martínez-Magaña JJ et al. (2023). Decoding the role of transcriptomic clocks in the human prefrontal cortex. medRxiv 2023.04.19.23288765; doi: https://doi.org/10.1101/2023.04.19.23288765
```

## Contact

For questions or contributions, contact:
**José Jaime Martínez-Magaña**
