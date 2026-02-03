# PRIME: scalable, robust inference of mechanistic cell states from multimodal single-cell counts via probability generating functions

## Introduction

In multimodal single-cell sequencing, cell-state heterogeneity is reflected not only in expression abundance but also in transcriptional dynamics and extrinsic variability. PRIME is a probability generating function (PGF)–based framework that infers mechanistic cell states from spliced/unspliced count data by (i) performing PGF-driven clustering and (ii) conducting Bayesian inference of gene- and cluster-specific kinetic parameters for downstream interpretation in parameter space. This repository (**Prime-analysis**) organizes the datasets and runnable notebooks used to reproduce the experiments and figures in the paper.

---

## Workflow

1. **dataset/** contains all datasets used in the manuscript, organized into three subfolders.
   - **synthetic/**: simulation datasets for **Fig. 2**.
     - **5-cluster setting**: includes files for `splice`, `unsplice`, `beta`, `celltype`, and ground-truth parameters `theta`.
     - **10-cluster setting**: includes files for `splice`, `unsplice`, `beta`, `celltype`, and ground-truth parameters `theta`.
     These datasets are used to evaluate performance and parameter recovery under controlled conditions.
   - **real_data/**: real datasets for **Fig. 3**.
     Real datasets used in **Fig. 3**. This folder contains files associated with five datasets:`cl3`, `cl5`, `Mopsc`, `humanskin`, and `mouselung`.For each dataset, we provide two types of inputs:
      - **Unfiltered `.loom` files** (e.g., `cl3_1137all.loom`, `cl5_1193all.loom`, `Human_skin_all.loom`, `mouse_all.loom`, `allen_b08_1948all.loom`):  
      These contain the complete count matrices and are used to estimate the **extrinsic-noise / cell-volume factor** \( \beta \).
      - **HVG-filtered count matrices in `.csv` format**:  
     These files contain **spliced** and **unspliced** count matrices restricted to **highly variable genes (HVGs)**, and are used as the primary inputs for PRIME clustering and downstream evaluation.
     Typical filenames include:
     - `*_hvg_spliced.csv`, `*_hvg_unspliced.csv` (HVG-filtered spliced/unspliced)
     - `spliced_expression_matrix_*.csv`, `unspliced_expression_matrix_*.csv` (HVG-filtered spliced/unspliced, dataset-specific naming)
     - `subclass_labels_*.csv` (ground-truth labels used for evaluation, e.g., ARI/NMI)

    In summary, **`.loom` files provide the full data required for estimating \( \beta \)**, whereas the **HVG-filtered `.csv` matrices are used for PRIME analysis** (clustering and metric reporting).
     - **mouselung**: input files for the mouse lung dataset.
   - **analysis/**: datasets for parameter-space marker analysis (**Figs. 5–6**).
     - **brca**: dataset corresponding to **Fig. 5** (BRCA analysis).
     - **ctrl** and **CTCL**: datasets corresponding to **Fig. 6**.
     - These datasets are used for marker gene analysis in the inferred parameter space, including volcano-plot style visualization of parameter-based markers.
     - **Note:** extremely large intermediate or posterior tables may be excluded from Git tracking due to file-size limits and should be generated locally by running the notebooks.

2. **run_prime/** contains runnable Jupyter notebooks (`.ipynb`) that reproduce the experiments.
   - **cl3.ipynb**: runs `prime_cluster` on the **cl3** dataset and reports clustering metrics (**ARI**, **NMI**).
   - **cl5.ipynb**: runs `prime_cluster` on the **cl5** dataset and reports (**ARI**, **NMI**).
   - **Mopsc.ipynb**: runs `prime_cluster` on **Mopsc** and reports (**ARI**, **NMI**).
   - **humanskin.ipynb**: runs `prime_cluster` on **humanskin** and reports (**ARI**, **NMI**).
   - **mouselung.ipynb**: runs `prime_cluster` on **mouselung** and reports (**ARI**, **NMI**).
   - **brca.ipynb**: reproduces the BRCA analysis workflow.
     - runs `prime_cluster` for clustering and compares discovered clusters against literature-reported markers;
     - runs `prime_infer_to_csv` for genome-wide Bayesian parameter inference;
     - performs parameter-space marker analysis and visualizes results (including marker volcano plots).

3. **(Dependency) Prime.jl** is required to execute the notebooks.
   - Install PRIME’s Julia package directly from GitHub:
     - `using Pkg`
     - `Pkg.develop(PackageSpec(url="https://github.com/Li-shiyue/Prime.jl.git"))`
     - `using Prime`
   - The notebooks in `run_prime/` assume the above installation is available in the Julia environment used by Jupyter.

---

## Repository Structure

```text
Prime-analysis/
  dataset/
    synthetic/      # Fig. 2
    real_data/      # Fig. 3
    analysis/       # Figs. 5–6
  run_prime/
    cl3.ipynb
    cl5.ipynb
    Mopsc.ipynb
    humanskin.ipynb
    mouselung.ipynb
    brca.ipynb
