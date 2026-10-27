# Structural Equation Modeling & Statistical Audit Pipeline (`R/`)

This directory contains `watershed_sem_pipeline.R`, the master statistical modeling script implementing all analyses reported in **"Microstructural Drivers of Lifespan Fluid Intelligence"**.

---

## 1. Execution Environment & Reproducibility

- **R Version**: `4.5.3`
- **Core Packages**: `lavaan` (`0.7-2`; Rosseel, 2012), `dplyr` (`1.2.1`), `parallel` (`mclapply`)
- **Parallelization & Determinism**: Uses `RNGkind("L'Ecuyer-CMRG")` seeded at `20261015L` across `19` CPU cores (`mc.cores = 19`) so that all `B = 5,000` parallel non-parametric bootstrap replicates are 100% deterministic and reproducible.
- **Pre-Filter Standardization Scope**: Because manifest variables are $z$-scored on $N = 547$ prior to excluding the 34 outliers/incomplete cases, the $N = 513$ analysis sample exhibits minor mean offsets ($M = -0.028\text{ to }+0.117$) and non-unit variances ($\text{Var} = 0.626\text{--}0.979$). All `lavaan` models therefore estimate full mean structures (`meanstructure = TRUE`, `fixed.x = FALSE`), extract completely standardized solutions (`standardizedSolution`), and use scale-aware standard-deviation weighting in all Wald equality tests.

---

## 2. Hierarchical Watershed SEM Specification ($df = 60$)

Across $p = 20$ manifest variables (4 Cattell indicators, 6 processing speed primitives, 10 bilateral white matter tracts), the forward watershed model specifies three structural tiers:

$$C_m = \nu_m + \lambda_m \text{FI} + \delta_m \quad (m = 1, \dots, 4)$$
$$\text{FI} = \sum_{j=1}^{6} \gamma_j \text{PS}_j + \zeta_{\text{FI}}$$
$$\text{PS}_j = \alpha_j + \sum_{k=1}^{10} \beta_{jk} \text{Tract}_k + \epsilon_j \quad (j = 1, \dots, 6)$$

- **Identification & Free Parameters**: Latent $\text{FI}$ is identified via marker-variable constraint ($\lambda_{C_1} = 1, \nu_{\text{FI}} = 0$). All $\binom{6}{2} = 15$ residual covariances among the processing speed primitives ($\epsilon_j \leftrightarrow \epsilon_{j'}$) and all $\binom{10}{2} = 45$ exogenous tract covariances ($\text{Tract}_k \leftrightarrow \text{Tract}_{k'}$, `fixed.x = FALSE`) are freely estimated alongside all 16 tract/speed variances and 20 manifest means/intercepts:
  $$\text{Sample Moments} = \frac{20 \times 21}{2} + 20 = 230, \qquad \text{Free Parameters} = 170, \qquad df = 230 - 170 = 60$$
- **Anatomy of the 60 Overidentifying Degrees of Freedom**: The 16-variable ($\text{WM}, \text{PS}$) covariance block is **fully saturated** with 136 covariance parameters ($16$ variances/residual variances + $\binom{16}{2} = 120$ off-diagonal covariances, regressions, and residual covariances) plus 16 means/intercepts. The 60 overidentifying degrees of freedom consist strictly of:
  1. **2 measurement-model tetrad constraints** among $C_1\text{--}C_4$
  2. **18 conditional-independence constraints** between the 6 PS primitives and 4 Cattell subtests ($C_m \perp \text{PS}_j \mid \text{FI}$)
  3. **40 conditional-independence constraints** between the 10 WM tracts and 4 Cattell subtests (10 direct $\text{WM} \rightarrow \text{FI}$ paths and 30 indicator-specific $\text{WM} \rightarrow C_m$ paths fixed to zero, enforcing $\text{WM} \perp \text{FI} \mid \text{PS}$ and $C_m \perp \text{WM} \mid \text{PS}$)
- **Estimator & Multiplicity Control**: Estimated via Robust Maximum Likelihood (`estimator = "MLR"`) with Yuan–Bentler scaled test statistics ($\chi^2_{\text{scaled}}$) and Huber–White robust standard errors. Within each model, the 60 upstream $\text{WM} \rightarrow \text{PS}$ $p$-values are adjusted using both the **Holm step-down FWER** ($p_{\text{Holm}}$) and **Benjamini–Hochberg FDR** ($p_{\text{BH}}$) procedures.

---

## 3. Mediator Ordering, Markov Equivalence, & Mediation Decomposition

For Models A and B, five structural audits are estimated on $N = 513$:

1. **Inverted Mediator Ordering ($\text{PS} \rightarrow \text{WM} \rightarrow \text{FI}$, $df = 56$, 174 free parameters)**: Regresses each of the 10 tracts on all 6 PS primitives (60 paths) and regresses $\text{FI}$ on the 10 tracts (10 paths, fixing direct $\text{PS} \rightarrow \text{FI}$ paths to zero). Because the 16-variable ($\text{WM}, \text{PS}$) block remains saturated, comparing $df = 60$ vs. $df = 56$ tests $\text{WM} \perp \text{FI} \mid \text{PS}$ against $\text{PS} \perp \text{FI} \mid \text{WM}$ (note: the $df = 56$ specification is Markov-equivalent to a common-cause topology $\text{PS} \leftarrow \text{WM} \rightarrow \text{FI}$ with $\text{PS} \leftrightarrow \text{PS}$).
2. **Reversed Lower Block ($\text{WM} \leftarrow \text{PS}$, $df = 60$)**: Reverses the 60 arrows between $\text{WM}$ and $\text{PS}$ while retaining $\text{FI} \sim \sum_j \text{PS}_j$ and fixing $\text{FI} \leftarrow 0 \cdot \text{Tract}_k$, verifying exact Markov equivalence ($\Delta\ell < 10^{-8}$, $\Delta\text{AIC} = 0.00$) of the saturated lower block.
3. **Partial-Mediation Model (`direct`, $df = 50$)**: Frees all 10 direct $\text{WM} \rightarrow \text{FI}$ paths and evaluates nested fit via the Satorra–Bentler scaled likelihood ratio test (`lavTestLRT`).
4. **Tract-Only Model (`wmonly`)**: Regresses $\text{FI}$ directly on the 10 tracts without PS variables to compute the proportion of WM-related variance in $\text{FI}$ mediated via PS:
   $$\text{Mediated Share} = \frac{R^2_{\text{wmonly}} - (R^2_{\text{direct}} - R^2_{\text{forward}})}{R^2_{\text{wmonly}}}$$
5. **Total Indirect Effects (`indirect`)**: Defines $\text{ind}_k := \sum_{j=1}^6 \beta_{jk}\gamma_j$ for each tract $k = 1, \dots, 10$, with delta-method standard errors on the completely standardized solution and 10-tract Holm correction.

---

## 4. Scale-Aware Wald Tests & Parallel Bootstrap Contrasts ($\text{ATR}$ vs. $\text{FMin}$)

- **Scale-Aware Wald Equality Tests (`lavTestWald`)**: To test $\beta_{j,\text{ATR}} = \beta_{j,\text{FMin}}$ within Models A and B while accounting for minor post-filter variance differences on $N = 513$, unstandardized regression weights are constrained proportional to their model-implied manifest standard deviations:
  $$\text{SD}_{\text{ATR}} \cdot b_{j,\text{ATR}} = \text{SD}_{\text{FMin}} \cdot b_{j,\text{FMin}} \quad (df = 1)$$
  with 6-outcome Holm adjustment ($p_{\text{Holm}}$) within each model.
- **Single-Pass Parallel Non-Parametric Bootstrap ($B = 5{,}000$)**: Resamples $N = 513$ participants with replacement (`100%` convergence across all 5,000 draws). Within each draw, unadjusted, `adj_age_sex`, and `adj_age_sex_vent` datasets are constructed and fitted across Models A and B to compute unadjusted marginal percentile $95\%$ CIs ($2.5\text{th}–97.5\text{th}$ empirical percentiles) across the 6 PS primitives for five focal contrasts (30 contrasts total):
  1. $\beta_{\text{ATR,A}} - \beta_{\text{FMin,A}}$
  2. $\beta_{\text{ATR,B}} - \beta_{\text{FMin,B}}$
  3. Cross-modality difference-in-differences: $(\beta_{\text{ATR,B}} - \beta_{\text{FMin,B}}) - (\beta_{\text{ATR,A}} - \beta_{\text{FMin,A}})$
  4. $\beta_{\text{FMin,B}} - \beta_{\text{FMin,A}}$
  5. $\beta_{\text{ATR,B}} - \beta_{\text{ATR,A}}$

---

## 5. Two-Stage Covariate Residualisation & Statistical Caveats

Prior to fitting the covariate-adjusted SEMs (`adj_age_sex` and `adj_age_sex_vent`), all 10 behavioral variables and 40 tract variables on $N = 513$ are residualised via ordinary least squares (QR decomposition) and re-standardized to unit variance:

- `adj_age_sex`: Partials out $\text{age}_z$, $\text{age}_z^2 = Z(\text{age}_z^2)$, and $\text{sex}_z$.
- `adj_age_sex_vent`: Additionally partials out $\text{vent}_z = Z(\ln[(V_{\text{Left Lat Vent}} + V_{\text{Right Lat Vent}})/\text{eTIV}])$.

**Important Psychometric Note (Manuscript Section 2.4.4 & Table 2)**: This two-stage OLS procedure partials out 3 to 4 covariates across all 20 manifest variables prior to SEM estimation, consuming 60 to 80 first-stage regression parameters not deducted from the second-stage SEM degrees of freedom ($df = 60$). This attenuates baseline off-diagonal correlations (lowering null-model $\chi^2$ and yielding $\text{TLI}_{\text{robust}} > 1.00$ whenever $\chi^2_{\text{scaled}} < df$) and alters the generalized variance $\det(\hat{\Sigma})$ upon unit-variance re-standardization. Consequently, joint AIC and BIC values are comparable only across models fitted to identical manifest variables and are not comparable between unadjusted and covariate-residualised regimes.
