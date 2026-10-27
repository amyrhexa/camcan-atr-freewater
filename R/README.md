# R Statistical Pipeline (`R/`)

This directory contains `watershed_sem_pipeline.R` (`wa.R`), a single-run, self-contained R script that executes all quality control, feature engineering, collinearity diagnostics, structural equation modeling (SEM), directionality/mediation audits, OLS covariate residualisation, parallel non-parametric bootstrapping, and outlier sensitivity analyses reported in the manuscript.

---

## Script Architecture (`watershed_sem_pipeline.R`)

The pipeline is organized into 10 sequential sections:

1. **Configuration, Path Resolution & Logging Engine:**
   - Forces single-threaded BLAS/OpenMP (`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`) so forked `parallel::mclapply` workers do not contend for CPU threads.
   - Auto-resolves the input dataset (`output/camcan_watershed_dataset.csv` or `camcan_watershed_dataset.csv`) and initializes deterministic parallel random number generation (`RNGkind("L'Ecuyer-CMRG")`, `seed = 20261015L`).
2. **Data Ingestion, Behavioural QC, Tract Aggregation & Cohort Profiling:**
   - Filters for valid behavioural task completion ($N = 547$), standardizes Cattell subtests (`C1`–`C4`), computes log-inverse speed and consistency primitives (`SRTspeed`, `SRTcons`, `CRTspeed`, `CRTcons`, `AVspeed`, `AVcons`), and aggregates bilateral JHU ROIs across Models A–D.
   - Applies the $\pm 4.0\text{ SD}$ outlier screen and completeness checks to derive the harmonized cohort ($N = 513$), exporting demographic summaries, age correlations, and Wilcoxon attrition audits (`excluded_vs_retained_age.csv`).
3. **Collinearity Diagnostics & Summary Covariance Deposits:**
   - Computes $10 \times 10$ inter-tract correlation matrices, Variance Inflation Factors ($\text{VIF}_k = (R^{-1})_{kk}$), and $20 \times 20$ manifest covariance matrices and mean vectors for Models A–D.
4. **Unified SEM Syntax Generator & Extraction Helpers:**
   - Programmatically constructs `lavaan` model syntax for 6 structural specifications: `forward` ($df = 60$), `direct` ($df = 50$), `indirect` ($df = 60$ with defined `:=` parameters), `reversed` ($df = 60$, $\text{WM} \leftarrow \text{PS}$), `inverted` ($df = 56$, $\text{PS} \rightarrow \text{WM} \rightarrow \text{FI}$), and `wmonly`.
   - Implements `extract_wald_tbl()` for scale-aware Wald tests (`lavTestWald`) testing standardized coefficient equality ($\text{SD}_{\text{ATR}} \cdot b_{j,\text{ATR}} = \text{SD}_{\text{FMin}} \cdot b_{j,\text{FMin}}$, $df = 1$) with 6-outcome Holm adjustment.
5. **Primary Watershed Models (A–D), Standalone CFA & Wald Tests:**
   - Estimates the standalone 1-factor Fluid Intelligence CFA ($df = 2$) and full forward watershed models A–D ($df = 60$) using Robust Maximum Likelihood (`estimator = "MLR"`, `missing = "fiml"`, `fixed.x = FALSE`).
   - Adjusts all 60 upstream $\text{WM} \rightarrow \text{PS}$ paths per model using Holm step-down FWER (`p_holm`) and Benjamini-Hochberg FDR (`p_bh`).
6. **Directionality, Hierarchy Inversion & Mediation Analyses:**
   - Evaluates the inverted mediator hierarchy ($\text{PS} \rightarrow \text{WM} \rightarrow \text{FI}$, $df = 56$), Markov equivalence of the lower-tier block ($\text{WM} \leftarrow \text{PS}$, $df = 60$), Satorra-Bentler scaled likelihood ratio tests (`lavTestLRT`) for freeing 10 direct $\text{WM} \rightarrow \text{FI}$ paths ($df = 50$), and delta-method total indirect effects ($\beta_{\text{ind},k} = \sum_{j=1}^6 \beta_{jk}\gamma_j$).
7. **Covariate-Residualised Watershed Models ($df = 60$):**
   - Performs vectorized QR-decomposition OLS residualisation (`qr.resid`) and unit-variance re-standardization across all 20 manifest variables under two regimes:
     - `adj_age_sex`: Controlling for standardized linear $\text{age}$, quadratic $Z(\text{age}^2)$, and $\text{sex}$.
     - `adj_age_sex_vent`: Additionally controlling for log-transformed FreeSurfer bilateral lateral ventricular volume normalized by $\text{eTIV}$:
       $$\text{vent}_z = Z\left(\ln\left(\frac{V_{\text{Left Lat Vent}} + V_{\text{Right Lat Vent}}}{\text{eTIV}}\right)\right)$$
8. **Single-Pass Parallel Bootstrap ($B = 5,000$):**
   - Uses `parallel::mclapply` and a fast ML point-estimate fitter (`se = "none"`, `test = "none"`) to resample $N = 513$ with replacement across all three specifications (`unadjusted`, `covariate_residualised_age_sex`, and `covariate_residualised_age_sex_vent`) simultaneously.
   - Extracts empirical percentile $95\%$ confidence intervals ($2.5\text{th}\text{--}97.5\text{th}$ percentiles) for 5 focal contrasts across all 6 speed outcomes:
     1. `ATR_minus_FMin_in_A`: $\beta_{\text{ATR},A} - \beta_{\text{FMin},A}$
     2. `ATR_minus_FMin_in_B`: $\beta_{\text{ATR},B} - \beta_{\text{FMin},B}$
     3. `change_in_contrast_B_vs_A` (Diff-in-Diff): $(\beta_{\text{ATR},B} - \beta_{\text{FMin},B}) - (\beta_{\text{ATR},A} - \beta_{\text{FMin},A})$
     4. `FMin_B_minus_A`: $\beta_{\text{FMin},B} - \beta_{\text{FMin},A}$
     5. `ATR_B_minus_A`: $\beta_{\text{ATR},B} - \beta_{\text{ATR},A}$
9. **Unified Outlier & Sample-Parity Sensitivity Module:**
   - Re-estimates Models A–D across three sensitivity regimes: `no_outlier_exclusion` ($N = 547$ Model A; $N = 538$ Models B–D), `modality_specific_4SD` ($N = 517\text{--}533$), and `common_noexcl` ($N = 538$).
10. **Session Info & Logging:**
    - Writes full R environment metadata to `sessionInfo.txt` and timestamped execution logs to `watershed_pipeline.log`.

---

## Usage & Environment Variables

| Environment Variable   | Default | Description                                                                                                             |
| :--------------------- | :------ | :---------------------------------------------------------------------------------------------------------------------- |
| `B_BOOT`               | `0`     | Number of non-parametric bootstrap replicates. Set `B_BOOT=5000` to reproduce the exact manuscript bootstrap intervals. |
| `OMP_NUM_THREADS`      | `1`     | OpenMP thread count (set to `1` for parallel worker safety).                                                            |
| `OPENBLAS_NUM_THREADS` | `1`     | OpenBLAS thread count (set to `1` for parallel worker safety).                                                          |
| `MKL_NUM_THREADS`      | `1`     | Intel MKL thread count (set to `1` for parallel worker safety).                                                         |

### Example Execution

```bash
# Full run with 5,000 parallel bootstrap resamples
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 B_BOOT=5000 Rscript R/watershed_sem_pipeline.R
```
