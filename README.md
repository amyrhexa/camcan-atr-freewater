# Microstructural Drivers of Lifespan Fluid Intelligence

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-4.5.3-276DC3.svg)](https://www.r-project.org/)
[![lavaan](https://img.shields.io/badge/lavaan-0.7--2-0f766e.svg)](https://lavaan.ugent.be/)
[![Python](https://img.shields.io/badge/Python-3.11-3776AB.svg)](https://www.python.org/)

Official analysis and visualization repository for the manuscript:

> **Microstructural Drivers of Lifespan Fluid Intelligence**
> Amir Hossein Mahmoudi$^{1,2}$, Zohreh Ganji$^{1,2}$, Arshia Alizadeh$^{1,2}$, Zeinab Mohammadpourmoghadam$^{1,2}$, Hoda Zare$^{1,3,*}$
> $^1$Department of Medical Physics, Faculty of Medicine, Mashhad University of Medical Sciences, Mashhad, Iran.
> $^2$Student Research Committee, Mashhad University of Medical Sciences, Mashhad, Iran.
> $^3$Medical Physics Research Center, Basic Sciences Research Institute, Mashhad University of Medical Sciences, Mashhad, Iran.
> $*$Corresponding author: `zareh@mums.ac.ir`

---

## Overview

Fluid intelligence declines steeply across the adult lifespan. The hierarchical **watershed model** conceptualizes fluid cognition as the downstream consequence of psychomotor processing speed endophenotypes governed by upstream cerebral white matter tracts. However, classical implementations relied on conventional single-tensor fractional anisotropy ($\text{FA}$), which conflates parenchymal axonal microstructure with isotropic extracellular free water ($\text{FW}$) driven by age-related ventriculomegaly.

In this study, we audit and extend the watershed framework across the adult lifespan ($N = 513$, age $18.6\text{--}88.0\text{ years}$) using multi-shell diffusion MRI ($b = 0, 1000, 2000\text{ s/mm}^2$) from the Cambridge Centre for Ageing and Neuroscience (**Cam-CAN**) cohort. Under strict sample parity ($N = 513$) and identical degrees of freedom ($df = 60$), we benchmark four biophysical diffusion representations:

- **Model A ($\text{DTI\_FA}$):** Conventional single-tensor fractional anisotropy ($b = 0, 1000\text{ s/mm}^2$).
- **Model B ($\text{FWE\_FA}_{\text{tissue}}$):** Bi-tensor free-water-eliminated tissue fractional anisotropy ($b = 0, 1000, 2000\text{ s/mm}^2$).
- **Model C ($\text{WMTI\_AWF}$):** White Matter Tract Integrity axonal water fraction ($b = 0, 1000, 2000\text{ s/mm}^2$).
- **Model D ($\text{FWE\_FW}$):** Isotropic extracellular free-water volume fraction ($b = 0, 1000, 2000\text{ s/mm}^2$).

### Key Empirical Findings

1. **Collinearity Reduction & Global Fit:** Bi-tensor $\text{FA}_{\text{tissue}}$ (Model B) achieves close approximate fit ($\text{RMSEA}_{\text{robust}} = 0.0325$, $\text{CFI}_{\text{robust}} = 0.9917$, $\text{SRMR} = 0.0298$) comparable to conventional $\text{FA}$ while substantially reducing inter-tract multicollinearity (maximum $\text{VIF}$ reduced from $3.226$ to $2.458$; inter-callosal correlation between Forceps Minor and Forceps Major reduced from $r = 0.656$ to $r = 0.350$).
2. **Proximate Mediator Ordering & Markov Equivalence:** Inverting the watershed hierarchy ($\text{PS} \rightarrow \text{WM} \rightarrow \text{FI}$, $df = 56$, enforcing $\text{PS} \perp \text{FI} \mid \text{WM}$) severely degrades model fit ($\Delta\text{AIC} = +160.26$ in Model B; $+118.48$ in Model A), confirming processing speed as the proximate mediator of fluid reasoning ($\text{WM} \perp \text{FI} \mid \text{PS}$). Reversing only the saturated 60-path lower block ($\text{WM} \leftarrow \text{PS}$, $df = 60$) is Markov-equivalent ($\Delta\text{AIC} = 0.00$).
3. **Shift in Anterior White Matter Specificity:** Eliminating extracellular free water resolves periventricular CSF partial-volume confounding in the callosal genu (Forceps Minor, $\text{FMin}$) and shifts upstream structural specificity to the **Anterior Thalamic Radiations ($\text{ATR}$)** (operationalized via bilateral anterior corona radiata ROIs in the JHU ICBM-DTI-81 atlas). In Model B, $\text{ATR} \rightarrow \text{CRT}_{\text{speed}}$ ($\beta = +0.356$, $p_{\text{Holm}} = 9.31 \times 10^{-8}$) and $\text{ATR} \rightarrow \text{CRT}_{\text{cons}}$ ($\beta = +0.277$, $p_{\text{Holm}} = 8.52 \times 10^{-4}$) are the sole upstream paths surviving 60-test family-wise error correction, and $\text{ATR}$ carries the sole Holm-significant total indirect effect on latent fluid intelligence ($\beta_{\text{ind}} = +0.269$, $p_{\text{Holm}} = 4.42 \times 10^{-8}$).
4. **Lifespan Covariance vs. Age-Invariant Variance:** OLS residualisation for $\text{age}$, $\text{age}^2$, $\text{sex}$, and FreeSurfer lateral ventricular volume normalized by intracranial volume ($\text{eTIV}$) preserves global model fit and downstream $\text{PS} \rightarrow \text{FI}$ coupling ($R^2 = 27.53\%$ vs. $57.77\%$ unadjusted) while attenuating upstream $\text{WM} \rightarrow \text{PS}$ regressions ($R^2 \le 5.43\%$), demonstrating that $\text{ATR}$ tissue anisotropy primarily mediates cross-sectional lifespan covariance.

---

## Repository Structure

```text
├── data/             # Input datasets, preprocessed SEM-ready data, and covariance deposits
│   └── README.md
├── R/                # End-to-end R Structural Equation Modeling (SEM) pipeline
│   ├── watershed_sem_pipeline.R
│   └── README.md
├── results/          # Generated outputs from the statistical and visualization pipelines
│   ├── figures       # 12 standalone vector publication figures (.pdf and .svg)
│   ├── tables/       # 28 exported statistical summary tables, fit indices, and logs (.csv)
│   └── README.md
├── scripts/          # Python publication figure generation workflow
│   ├── generate_publication_figures.py
│   └── README.md
└── README.md         # Project overview and reproduction guide
```

---

## System Requirements & Dependencies

### R Environment (`R/watershed_sem_pipeline.R`)

- **R** $\ge 4.5.0$ (tested on **v4.5.3**)
- Required packages:
  - `lavaan` (tested on **v0.7-2**)
  - `dplyr`
  - `parallel` (base R)

Install required R packages via:

```R
install.packages(c("dplyr", "lavaan"))
```

### Python Environment (`scripts/generate_publication_figures.py`)

- **Python** $\ge 3.10$ (tested on **v3.11**)
- Required packages:
  - `numpy`
  - `pandas`
  - `matplotlib`
  - `seaborn`

Install required Python packages via:

```bash
pip install numpy pandas matplotlib seaborn
```

---

## End-to-End Reproduction Guide

### Step 1: Prepare Input Data

Place the merged Cam-CAN Stage 2 dataset (`camcan_watershed_dataset.csv`) in the expected input directory (see `data/README.md` for required schema and access details).

_Note:_ `R/watershed_sem_pipeline.R` automatically checks `output/camcan_watershed_dataset.csv` or `./camcan_watershed_dataset.csv`. If your raw file is stored in `data/camcan_watershed_dataset.csv`, create a symbolic link or copy before running:

```bash
mkdir -p output
ln -sf "$(pwd)/data/camcan_watershed_dataset.csv" output/camcan_watershed_dataset.csv
```

### Step 2: Run the Watershed SEM Pipeline (R)

Execute `R/watershed_sem_pipeline.R` (`wa.R`). To reproduce the full manuscript analysis including the $B = 5,000$ single-pass parallel non-parametric bootstrap (with L'Ecuyer-CMRG seed `20261015L` and single-threaded BLAS per worker to prevent thread contention):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 B_BOOT=5000 Rscript R/watershed_sem_pipeline.R
```

_Quick Run (Point Estimates, Wald Tests, Directionality, and Sensitivity without Bootstrap):_

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 B_BOOT=0 Rscript R/watershed_sem_pipeline.R
```

Sync the exported CSV tables and logs into `results/tables/`:

```bash
mkdir -p results/tables
cp output/*.csv output/*.txt output/*.log results/tables/
```

### Step 3: Generate Publication Figures (Python)

Run `scripts/generate_publication_figures.py` pointing `--data_dir` to `results/tables` and `--output_dir` to `results/figures`:

```bash
python3 scripts/generate_publication_figures.py \
  --data_dir results/tables \
  --output_dir results/figures
```

This script generates all 12 standalone vector figures in both `.pdf` and `.svg` formats inside `results/figures/`.

---

## Data & Code Availability

- **Data Availability:** Raw and scored behavioural and neuroimaging data from the Cambridge Centre for Ageing and Neuroscience (Cam-CAN) are available to researchers upon application via the [Cam-CAN Data Access Portal](https://www.cam-can.com) (Shafto et al., _BMC Neurol._, 2014). Summary covariance matrices and variable means for Models A–D ($N = 513$) are exported directly by the R pipeline (`covariance_matrix_A.csv`–`D.csv` and `means_A.csv`–`D.csv`), allowing full SEM replication in `lavaan` without subject-level data.
- **Code Availability:** All R and Python workflows in this repository are released under the **MIT License** at [https://github.com/amyrhexa/camcan-atr-freewater](https://github.com/amyrhexa/camcan-atr-freewater).

---

## Citation

If you use this codebase, pipeline architecture, or summary covariance matrices in your research, please cite:

```bibtex
@article{mahmoudi2026microstructural,
  title   = {Microstructural Drivers of Lifespan Fluid Intelligence},
  author  = {Mahmoudi, Amir Hossein and Ganji, Zohreh and Alizadeh, Arshia and Mohammadpourmoghadam, Zeinab and Zare, Hoda},
  year    = {2026},
  note    = {Submitted manuscript},
  url     = {https://github.com/amyrhexa/camcan-atr-freewater}
}
```
