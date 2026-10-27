# Microstructural Drivers of Lifespan Fluid Intelligence

[![Preprint: NeuroImage](<https://img.shields.io/badge/Preprint-NeuroImage_(2026)-0056b3.svg>)](https://github.com/amyrhexa/camcan-atr-freewater)
[![R](https://img.shields.io/badge/R-4.5.3-276DC3.svg?logo=r&logoColor=white)](https://www.r-project.org/)
[![lavaan](https://img.shields.io/badge/lavaan-0.7--2-4B8BBE.svg)](https://lavaan.ugent.be/)
[![Python](https://img.shields.io/badge/Python-3.11-3776AB.svg?logo=python&logoColor=white)](https://www.python.org/)
[![DIPY](https://img.shields.io/badge/DIPY-1.12.1-FF6F00.svg)](https://dipy.org/)

Official statistical modeling, biophysical auditing, and visualization repository for the manuscript:

# **Microstructural Drivers of Lifespan Fluid Intelligence**

> **Amir Hossein Mahmoudi**<sup>a,b</sup>, **Zohreh Ganji**<sup>a,b</sup>, **Arshia Alizadeh**<sup>a,b</sup>, **Zeinab Mohammadpourmoghadam**<sup>a,b</sup>, **Hoda Zare**<sup>a,c,\*</sup>

> <sup>a</sup> _Department of Medical Physics, Faculty of Medicine, Mashhad University of Medical Sciences, Mashhad, Iran_

> <sup>b</sup> _Student Research Committee, Mashhad University of Medical Sciences, Mashhad, Iran_

> <sup>c</sup> _Medical Physics Research Center, Basic Sciences Research Institute, Mashhad University of Medical Sciences, Mashhad, Iran_

> \*Corresponding author: `zareh@mums.ac.ir` (Hoda Zare)

> \*Preprint submitted to NeuroImage (November 1, 2026)\*

---

## Overview

The hierarchical **watershed model** of cognitive aging (Cannon & Keller, 2006; Kievit et al., 2016) posits a three-tier many-to-one cascade across the adult lifespan: multiple cerebral white matter (WM) tracts govern distinct psychomotor processing speed (PS) primitives, which in turn funnel variance to downstream latent fluid intelligence (FI). However, classical watershed implementations relied exclusively on single-tensor fractional anisotropy ($\text{FA}$), which is confounded at standard 2 mm isotropic resolution by age-related extracellular isotropic free water ($\text{FW}$) from ex-vacuo ventriculomegaly, sulcal widening, and perivascular space dilation.

Because the genu of the corpus callosum (**Forceps Minor**, $\text{FMin}$) forms the anterior wall and roof of the frontal horns of the lateral ventricles, age-related ventricular enlargement exposes $\text{FMin}$ to pronounced cerebrospinal fluid (CSF) partial volume effects. We audited and extended the watershed structural equation modeling (SEM) framework across the adult lifespan ($N = 513$ harmonized complete cases; age $18.6\text{--}88.0\text{ years}$, $253\text{ females}$, $260\text{ males}$) using multi-shell diffusion MRI ($b = 0, 1000, 2000\text{ s/mm}^2$) from the Cambridge Centre for Ageing and Neuroscience (**Cam-CAN** Stage 2 cohort).

---

## Core Empirical Findings

### 1. Multimodal Biophysical Benchmark & Collinearity Reduction ($N = 513$, $df = 60$)

Four parallel watershed SEMs were estimated on the identical $N = 513$ cohort using Robust Maximum Likelihood (`MLR`, Yuan–Bentler scaled $\chi^2_{\text{scaled}}$, Huber–White robust standard errors):

- **Model B (Bi-Tensor Free-Water-Eliminated Tissue FA, $\text{FA}_{\text{tissue}}$)**: $\chi^2_{\text{scaled}} = 92.270$ ($p = 0.0047$), $\text{CFI}_{\text{robust}} = 0.9917$, $\text{TLI}_{\text{robust}} = 0.9737$, $\text{RMSEA}_{\text{robust}} = 0.0325$ ($90\%\text{ CI }[0.0183, 0.0452]$), $\text{SRMR} = 0.0298$.
- **Model D (Extracellular Free-Water Fraction, $\text{FW}$)**: $\chi^2_{\text{scaled}} = 91.677$ ($p = 0.0053$), $\text{CFI}_{\text{robust}} = 0.9934$, $\text{TLI}_{\text{robust}} = 0.9791$, $\text{RMSEA}_{\text{robust}} = 0.0332$ ($90\%\text{ CI }[0.0176, 0.0468]$), $\text{SRMR} = 0.0313$.
- **Model A (Conventional Single-Shell FA, $b = 1000\text{ s/mm}^2$)**: $\chi^2_{\text{scaled}} = 97.561$ ($p = 0.0016$), $\text{CFI}_{\text{robust}} = 0.9913$, $\text{TLI}_{\text{robust}} = 0.9724$, $\text{RMSEA}_{\text{robust}} = 0.0354$ ($90\%\text{ CI }[0.0219, 0.0479]$), $\text{SRMR} = 0.0333$.
- **Model C (WMTI Axonal Water Fraction, $\text{AWF}$)**: $\chi^2_{\text{scaled}} = 108.764$ ($p = 0.0001$), $\text{CFI}_{\text{robust}} = 0.9903$, $\text{TLI}_{\text{robust}} = 0.9694$, $\text{RMSEA}_{\text{robust}} = 0.0401$ ($90\%\text{ CI }[0.0279, 0.0520]$), $\text{SRMR} = 0.0372$.

_Methodological Note_: Because the $90\%\text{ RMSEA}_{\text{robust}}$ confidence intervals overlap substantially across Models A–D ($\Delta\text{RMSEA}_{\text{robust}} = 0.0029$ between Models B and A) and each model is fitted to distinct modality-specific observed white matter variables, small descriptive differences in global fit indices are not inferentially significant. Rather, the primary methodological advantage of Model B ($\text{FA}_{\text{tissue}}$) is its biophysical separation of parenchymal tissue anisotropy from extracellular free water and its substantially lower inter-tract multicollinearity: maximum tract Variance Inflation Factor ($\text{VIF}$) dropped from **4.065** ($\text{AWF}$), **3.829** ($\text{FW}$), and **3.226** (conventional $\text{FA}$) to **2.458** in $\text{FA}_{\text{tissue}}$, with the inter-callosal correlation between $\text{FMin}$ and $\text{FMaj}$ declining from $r = 0.656$ (Model A) to $r = 0.350$ (Model B).

### 2. Mediator Ordering, Markov Equivalence, and Variance Decomposition

- **Rejection of Inverted Mediator Ordering ($\text{PS} \rightarrow \text{WM} \rightarrow \text{FI}$, $df = 56$)**: Enforcing conditional independence $\text{PS} \perp \text{FI} \mid \text{WM}$ severely degraded fit in both Model B ($\chi^2_{\text{scaled}} = 241.169$, $\text{CFI}_{\text{robust}} = 0.9520$, $\text{RMSEA}_{\text{robust}} = 0.0810$, $\Delta\text{AIC} = +160.26$) and Model A ($\chi^2_{\text{scaled}} = 205.519$, $\Delta\text{AIC} = +118.48$), confirming processing speed as the proximate mediator of fluid intelligence ($\text{WM} \perp \text{FI} \mid \text{PS}$).
- **Lower-Block Markov Equivalence ($df = 60$)**: Because the 16-variable ($\text{WM}, \text{PS}$) covariance block is fully saturated (136 covariance parameters + 16 means), reversing only the 60 arrows between $\text{WM}$ and $\text{PS}$ ($\text{WM} \leftarrow \text{PS}$ while retaining $\text{PS} \rightarrow \text{FI}$ and $\text{FI} \leftarrow 0 \cdot \text{WM}$, $df = 60$) is Markov-equivalent to the forward model ($\Delta\ell < 10^{-8}$, $\Delta\text{AIC} = 0.00$). Orienting $\text{WM} \rightarrow \text{PS}$ relies on neurobiological and developmental priors.
- **Mediation Proportion & Absence of Direct $\text{FMin} \rightarrow \text{FI}$ Path**: The six processing speed primitives explained **$57.77\%$** of variance in latent $\text{FI}$ ($\text{CRT}_{\text{speed}}\;\beta = +0.564$, $\text{AV}_{\text{cons}}\;\beta = +0.206$, $\text{CRT}_{\text{cons}}\;\beta = +0.117$, $\text{SRT}_{\text{cons}}\;\beta = +0.093$, suppressor $\text{AV}_{\text{speed}}\;\beta = -0.252$, non-significant $\text{SRT}_{\text{speed}}\;\beta = -0.034$). Freeing all 10 direct $\text{WM} \rightarrow \text{FI}$ paths ($df = 50$) increased $R^2(\text{FI})$ in Model B from $57.77\%$ to $62.19\%$ ($+4.42$ percentage points vs. $R^2_{\text{wmonly}} = 36.10\%$), indicating that **$87.8\%$** of white-matter-related variance in $\text{FI}$ is shared with processing speed. Crucially, the direct $\text{FMin} \rightarrow \text{FI}$ path was non-significant in both Model A ($\beta = +0.120$, $p = 0.0517$) and Model B ($\beta = +0.085$, $p = 0.0736$), and BIC favored the parsimonious full-mediation model ($\Delta\text{BIC} = +24.13$ in Model B; $+15.51$ in Model A).

### 3. Shift in Anterior Specificity: Forceps Minor ($\text{FMin}$) vs. Anterior Thalamic Radiations ($\text{ATR}$)

- **Conventional FA (Model A)**: $\text{FMin}$ and $\text{ATR}$ were virtually indistinguishable for $\text{CRT}_{\text{speed}}$ ($\beta_{\text{ATR}} = +0.182$ vs. $\beta_{\text{FMin}} = +0.181$; scale-aware Wald $p = 0.9906$), $\text{FMin}$ had the largest point estimate for $\text{CRT}_{\text{cons}}$ ($\beta_{\text{FMin}} = +0.243$ vs. $\beta_{\text{ATR}} = +0.127$), and no upstream path survived 60-test Holm FWER correction ($p_{\text{Holm}} \ge 0.0562$).
- **Periventricular Free-Water Mechanism (Model D)**: Extracellular free water in $\text{FMin}$ negatively predicted $\text{CRT}_{\text{cons}}$ ($\beta = -0.213$, $p = 0.0090$) and marginally predicted $\text{CRT}_{\text{speed}}$ ($\beta = -0.121$, $p = 0.0563$), whereas $\text{ATR}$ free water had near-zero effects on both choice primitives ($\beta = -0.022$ and $-0.024$, $p > 0.73$).
- **Free-Water-Eliminated Tissue FA (Model B)**: Removing extracellular free water elevated **$\text{ATR}$** (operationalized in the JHU ICBM-DTI-81 atlas via the bilateral anterior corona radiata ROIs) to the sole white matter tract surviving 60-path Holm FWER and Benjamini–Hochberg FDR corrections against zero:
  - $\text{ATR} \rightarrow \text{CRT}_{\text{speed}}$: $\beta = +0.356$ ($z = 6.04$, $p = 1.55 \times 10^{-9}$, $p_{\text{Holm}} = 9.31 \times 10^{-8}$, $p_{\text{BH}} = 9.31 \times 10^{-8}$) vs. $\text{FMin}$ $\beta = +0.115$ ($p = 0.0251$, $p_{\text{Holm}} = 1.00$).
  - $\text{ATR} \rightarrow \text{CRT}_{\text{cons}}$: $\beta = +0.277$ ($z = 4.34$, $p = 1.44 \times 10^{-5}$, $p_{\text{Holm}} = 8.52 \times 10^{-4}$, $p_{\text{BH}} = 4.33 \times 10^{-4}$) vs. $\text{FMin}$ $\beta = +0.134$ ($p = 0.0356$, $p_{\text{Holm}} = 1.00$).
  - **Total Indirect Effect on FI**: Only $\text{ATR}$ survived 10-tract Holm correction in Model B ($\beta_{\text{ind}} = +0.269$, $z = 5.87$, $p = 4.42 \times 10^{-9}$, $p_{\text{Holm}} = 4.42 \times 10^{-8}$), $2.49\times$ larger than $\text{FMin}$ ($\beta_{\text{ind}} = +0.108$, $p = 0.0130$, $p_{\text{Holm}} = 0.1167$).
- **Formal Equality & Cross-Modality Bootstrap Contrasts ($B = 5{,}000$)**:
  - _Cross-Modality Difference-in-Differences_: Paired bootstrap contrasts $\Delta\Delta\beta = (\beta_{\text{ATR,B}} - \beta_{\text{FMin,B}}) - (\beta_{\text{ATR,A}} - \beta_{\text{FMin,A}})$ showed positive shifts toward $\text{ATR}$ relative to $\text{FMin}$ whose unadjusted marginal $95\%$ percentile CIs excluded zero for both $\text{CRT}_{\text{speed}}$ ($+0.2404\;[+0.0996, +0.3904]$) and $\text{CRT}_{\text{cons}}$ ($+0.2579\;[+0.1117, +0.4071]$).
  - _Within-Model B Equality Caution_: Direct within-Model B scale-aware Wald equality tests ($\beta_{\text{ATR,B}} = \beta_{\text{FMin,B}}$) were nominally significant for $\text{CRT}_{\text{speed}}$ ($\Delta\beta = +0.2416$, $\chi^2(1) = 5.943$, $p = 0.0148$) but did not survive 6-outcome Holm correction ($p_{\text{Holm}} = 0.0887$), and were non-significant for $\text{CRT}_{\text{cons}}$ ($\Delta\beta = +0.1421$, $\chi^2(1) = 1.472$, $p = 0.2251$) and the remaining four primitives ($p \ge 0.4791$), reflecting shared variance between $\text{ATR}$ and $\text{FMin}$ ($r = 0.701$).

### 4. Lifespan Covariance vs. Age- and Ventriculomegaly-Residualised Regimes

Controlling for linear $\text{age}$, quadratic $\text{age}^2$, $\text{sex}$, and FreeSurfer normalized bilateral lateral ventricular volume ($\ln[(\text{Left}+\text{Right Lat Vent})/\text{eTIV}]$) preserved close global fit ($\chi^2_{\text{scaled}} = 54.660\text{--}64.123$, $df = 60$, $p = 0.3341\text{--}0.6705$) and substantial $\text{PS} \rightarrow \text{FI}$ coupling ($R^2 = 27.51\%\text{--}27.53\%$), but attenuated upstream $\text{WM} \rightarrow \text{PS}$ variance ($R^2 = 5.43\%$ for $\text{CRT}_{\text{speed}}$; $\text{ATR} \rightarrow \text{CRT}_{\text{speed}}\;\beta = +0.045$, $p = 0.4391$, though the paired cross-modality shift $\Delta\beta_{\text{ATR,B}-\text{A}} = +0.0576\;[+0.0025, +0.1137]$ remained positive). Thus, $\text{ATR}$ tissue anisotropy primarily mediates cross-sectional lifespan aging covariance, whereas age-invariant effects on simple/bimodal tempo are modest and localized to the Cingulate Cingulum ($\text{CING}$).

---

## Repository Structure

```text
camcan-atr-freewater/
├── README.md                                 # Project overview, empirical summary, and reproduction guide
├── CITATION.cff                              # Machine-readable citation metadata
├── LICENSE                                   # MIT License
├── data/
│   └── README.md                             # Cam-CAN Stage 2 cohort harmonization, QC rules, and variable codebook
├── R/
│   ├── README.md                             # Statistical methodology & lavaan SEM specification guide
│   └── watershed_sem_pipeline.R              # Master R script: cohort audit, SEMs, equivalence tests, 5,000-rep bootstrap
├── scripts/
│   ├── README.md                             # Figure generation documentation & manuscript figure mapping
│   └── generate_publication_figures.py       # Python script generating all 12 vector PDF/SVG manuscript figures
└── results/
    ├── README.md                             # Complete index mapping 38 CSV tables & 12 figures to the manuscript
    ├── figures/                              # 12 publication figures in vector (.pdf) and web (.svg) formats
    ├── tables/                               # 38 structured CSV tables (fits, parameters, VIFs, Wald & bootstrap CIs)
    ├── watershed_pipeline.log                # Execution log from watershed_sem_pipeline.R
    └── sessionInfo.txt                       # Exact R 4.5.3 environment & package state
```

---

## Reproducing the Analysis & Figures

### Software Requirements

- **Neuroimaging Preprocessing Stack**: MRtrix3 (`3.0.8`), FSL `eddy` (`6.0.7.22`), FreeSurfer (`8.2.0`), DIPY (`1.12.1`) in Python (`3.11`).
- **Statistical Modeling Environment**: R (`4.5.3`) with `lavaan` (`0.7-2`), `dplyr` (`1.2.1`), and base `parallel` (`mclapply`, 19 cores, RNG `L'Ecuyer-CMRG` seeded at `20261015L`).
- **Visualization Environment**: Python (`3.11+`) with `numpy`, `pandas`, `matplotlib`, and `seaborn`.

### Step 1: Run the Structural Equation Modeling Pipeline

Once the Cam-CAN Stage 2 source tables are placed in `data/` (see [`data/README.md`](data/README.md)), execute the master R pipeline from the repository root:

```bash
Rscript R/watershed_sem_pipeline.R
```

This script executes the 4-stage cohort harmonization ($N = 653 \rightarrow 550 \rightarrow 547 \rightarrow 513$), runs the Wilcoxon attrition audit, fits Models A–D across unadjusted, age/sex-residualised, and age/sex/ventricle-residualised regimes, runs the mediator-ordering and Markov-equivalence models, executes the 5,000-replicate parallel non-parametric bootstrap, evaluates the three outlier sensitivity regimes ($N = 517\text{--}547$), and writes all 38 CSV tables to `results/tables/`.

### Step 2: Regenerate All 12 Manuscript Figures

```bash
python3 scripts/generate_publication_figures.py
```

Reads the summary CSVs in `results/tables/` and outputs all 12 figures (`fig01`–`fig12` in `.pdf` and `.svg`) to `results/figures/`.

---

## Citation

If you use this pipeline, statistical audit framework, or codebase in your research, please cite:

```bibtex
@article{mahmoudi2026microstructural,
  title   = {Microstructural Drivers of Lifespan Fluid Intelligence},
  author  = {Mahmoudi, AmirHossein and Ganji, Zohreh and Alizadeh, Arshia and Mohammadpourmoghadam, Zeinab and Zare, Hoda},
  journal = {NeuroImage},
  year    = {2026},
  note    = {Preprint submitted November 1, 2026},
  url     = {https://github.com/amyrhexa/camcan-atr-freewater}
}
```

## License

Code in this repository is released under the MIT License. Raw and derived Cam-CAN data are subject to the data usage agreements of the Cambridge Centre for Ageing and Neuroscience (`https://www.cam-can.com`).
