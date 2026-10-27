# Data Access, Cohort Harmonization, and Variable Codebook

This directory houses the input data tables and harmonized analysis matrix (`camcan_watershed_dataset.csv`) for **"Microstructural Drivers of Lifespan Fluid Intelligence"** ($N = 513$).

---

## 1. Data Access & Ethics Statement

Raw neuroimaging, demographic, and behavioral data were obtained from the **Stage 2** cross-sectional adult lifespan repository of the **Cambridge Centre for Ageing and Neuroscience (Cam-CAN)** (Shafto et al., 2014).

- **Data Portal**: Researchers can request access via the official Cam-CAN portal at [https://www.cam-can.com](https://www.cam-can.com).
- **Ethics Approval**: Granted by the Cambridgeshire 2 Research Ethics Committee (reference: `10/H0308/50`). All participants provided written informed consent.
- **Inclusion Screening**: Community-dwelling adults aged 18–88 years screened to exclude Mini-Mental State Examination ($\text{MMSE}$) scores $< 24$, major neurological or psychiatric conditions, uncorrected sensory deficits, non-fluent English, or MRI contraindications.

Because Cam-CAN data usage agreements prohibit redistributing raw participant-level repository dumps without an approved DUA, place the licensed Cam-CAN Stage 2 CSV tables in this directory before executing `R/watershed_sem_pipeline.R`.

---

## 2. Four-Stage Cohort Derivation & Attrition Audit (Manuscript Figure 1)

Cohort harmonization follows a deterministic four-stage pipeline implemented in `R/watershed_sem_pipeline.R`:

1. **Stage 1 — Cross-Modal Sequential Inner Merge ($N = 653 \rightarrow N = 550$)**:
   - Core demographics (`participants.csv`, $N = 653$) and extended demography ($N = 2{,}681 \rightarrow 653$ retained)
   - Cattell Culture Fair Intelligence Scale 2 Form A ($N = 660 \rightarrow 634$ retained)
   - Simple Reaction Time ($N = 2{,}374 \rightarrow 588$ retained)
   - Choice Reaction Time ($N = 2{,}327 \rightarrow 577$ retained)
   - MRI behavioral Audio-Visual Cued Reaction Time ($N = 638 \rightarrow 564$ retained)
   - Multi-shell JHU ICBM-DTI-81 white matter ROI tables ($N = 635 \rightarrow 550$ retained)
   - FreeSurfer `aseg` morphometry (`global_summary.csv` and `aseg_volume.csv`, $N = 648 \rightarrow 550$ retained)
   - _Result_: Post-merge dataset of **$N = 550$ participants and 561 features**. An automated range check flags extreme unscaled numerical overflow in unused WMTI tortuosity columns ($\max = 4.13 \times 10^8$) and excludes them from all analyses.
2. **Stage 2 — Behavioral Quality Control ($N = 550 \rightarrow N = 547$)**:
   - Requires binary task-inclusion flags equal to `1` across all three reaction time paradigms: `rt_simple_Include_participant == 1`, `rt_choice_Include_participant == 1`, and `mri_behav_include_participant == 1`.
   - Excludes 3 participants, retaining **$N = 547$ behaviorally valid participants** ($99.45\%$ retention; $0$ missing behavioral observations).
3. **Stage 3 — Standardization, 4-SD Outlier Screening, & Multi-Compartment Completeness ($N = 547 \rightarrow N = 513$)**:
   - Across the $N = 547$ behaviorally valid cohort, all 10 behavioral indicators and 40 bilateral white matter tract variables ($10\text{ tracts} \times 4\text{ modalities}$) are standardized to $z$-scores ($M = 0.000, \text{SD} = 1.000$).
   - Following Kievit et al. (2016), participants with $|z| > 4.0$ on any manifest behavioral or tract indicator are flagged as univariate outliers: **25 unique participants** ($4.57\%$: 4 behavioral outliers and 22 tract outliers, with 1 overlapping both).
   - **9 participants** ($1.64\%$) exhibit missing values across the multi-compartment diffusion maps ($\text{FA}_{\text{tissue}}, \text{AWF}, \text{FW}$) due to non-convergence/boundary exclusion, whereas single-tensor $\text{FA}$ has $0$ missing cases.
4. **Stage 4 — Final Harmonized Complete-Case SEM Cohort ($N = 513$)**:
   - Excluding the 25 outliers and 9 incomplete multi-compartment dMRI cases yields **$N = 513$ participants** ($253\text{ females }[49.3\%], 260\text{ males }[50.7\%]$; age $M = 53.75\text{ years}, \text{SD} = 17.95, \text{range} = 18.6\text{--}88.0\text{ years}, \text{median} = 52.90, \text{IQR} = [38.90, 68.80]$).
   - **Wilcoxon Rank-Sum Attrition Audit**: The 34 excluded participants are significantly older than the 513 retained participants ($M = 70.55\text{ vs. }53.75\text{ years}, p = 1.19 \times 10^{-7}$), driven by both the 22 tract outliers ($M = 70.45\text{ years}, p = 2.31 \times 10^{-5}$) and the 9 incomplete multi-compartment dMRI cases ($M = 71.90\text{ years}, \text{range }69.6\text{--}75.5\text{ years}, p = 0.0021$). All models are therefore verified across three sensitivity regimes retaining the outliers ($N = 517\text{--}547$).

---

## 3. Neuroimaging Preprocessing & Biophysical Models

- **Acquisition**: 3T Siemens TIM Trio (32-channel head coil), twice-refocused spin-echo EPI dMRI ($\text{TR} = 9100\text{ ms}$, $\text{TE} = 104\text{ ms}$, $2\text{ mm}$ isotropic voxels, 30 directions at $b = 1000\text{ s/mm}^2$, 30 directions at $b = 2000\text{ s/mm}^2$, 3 interspersed $b = 0\text{ s/mm}^2$ volumes) + $1\text{ mm}$ isotropic T1w MPRAGE and T2w SPACE.
- **Preprocessing**: Marchenko–Pastur PCA thermal denoising and Gibbs unringing in **MRtrix3 (`3.0.8`)**, multi-shell eddy-current/motion correction with slice-to-volume outlier replacement and dynamic $b$-vector rotation in **FSL `eddy`**, ANTs N4 bias-field correction, and boundary-based registration (BBR) + `aseg` ventricular/eTIV morphometry in **FreeSurfer (`8.2.0`)**.
- **Voxelwise Biophysical Fitting in DIPY (`1.12.1`, Python `3.11`)**:
  1. **Model A ($\text{FA}$)**: Single-tensor DTI fitted to the $b = 0, 1000\text{ s/mm}^2$ shell for direct methodological parity with Kievit et al. (2016).
  2. **Model B ($\text{FA}_{\text{tissue}}$) & Model D ($\text{FW}$)**: Fitted across all shells ($b = 0, 1000, 2000\text{ s/mm}^2$) via the bi-tensor Free-Water Elimination model:
     $$S_k = S_0 \left[(1 - f_{\text{iso}})\exp(-b_k \mathbf{g}_k^\top \mathbf{D}_{\text{tissue}} \mathbf{g}_k) + f_{\text{iso}}\exp(-b_k d_{\text{water}})\right], \quad d_{\text{water}} = 3.0 \times 10^{-3}\text{ mm}^2/\text{s}$$
  3. **Model C ($\text{AWF}$)**: Derived across $b = 0, 1000, 2000\text{ s/mm}^2$ via the White Matter Tract Integrity (WMTI) kurtosis model:
     $$\text{AWF} = \frac{K_{\max}}{K_{\max} + 3}$$

---

## 4. Manifest Variable Codebook (Manuscript Table 1)

### Tier 1: Fluid Intelligence Indicators (Cattell Culture Fair Scale 2, Form A)

Standardized to $z$-scores ($M = 0, \text{SD} = 1$) on $N = 547$:

- `C1` — Series Completion ($r_{\text{age}} = -0.545$)
- `C2` — Classification ($r_{\text{age}} = -0.499$)
- `C3` — Matrices ($r_{\text{age}} = -0.584$)
- `C4` — Conditions / Topology ($r_{\text{age}} = -0.468$)

### Tier 2: Psychomotor Processing Speed Primitives

In raw Cam-CAN behavioral tables, SRT and CRT latencies are recorded in seconds, whereas AV trial means (`mri_behav_mRT`) and SDs (`mri_behav_stdRT`) are stored in milliseconds. Values $> 10$ are divided by $1{,}000$ to convert uniformly to seconds (all $N = 547$ participants have AV mean and SD $> 10\text{ ms}$, avoiding any threshold discontinuity). All latency means and intra-individual SDs are inverted, natural-log-transformed, and $z$-scored on $N = 547$ so higher values indicate **faster speed** and **greater consistency**:
$$\text{Speed}_j = Z\!\left(\ln\!\left(\frac{1}{\text{Mean RT}_j}\right)\right), \qquad \text{Consistency}_j = Z\!\left(\ln\!\left(\frac{1}{\text{SD RT}_j}\right)\right)$$

- `SRT_speed`, `SRT_cons` — Simple Reaction Time speed ($r_{\text{age}} = -0.367$) and consistency ($r_{\text{age}} = -0.330$)
- `CRT_speed`, `CRT_cons` — 4-Choice Reaction Time speed ($r_{\text{age}} = -0.657$) and consistency ($r_{\text{age}} = -0.558$)
- `AV_speed`, `AV_cons` — Unhurried Audio-Visual Cued speed ($r_{\text{age}} = +0.027$) and consistency ($r_{\text{age}} = -0.248$)

### Tier 3: 10 Bilateral JHU ICBM-DTI-81 White Matter Tract Systems

Aggregated across each modality ($\text{FA}$, $\text{FA}_{\text{tissue}}$, $\text{AWF}$, $\text{FW}$) and $z$-scored on $N = 547$:

1. `FMin` — **Forceps Minor**: Genu of the corpus callosum ($\text{FA}_{\text{tissue}}\;r_{\text{age}} = -0.584$; $\text{VIF}_{B} = 2.241$ vs. $\text{VIF}_{A} = 3.226$)
2. `FMaj` — **Forceps Major**: Splenium of the corpus callosum ($\text{FA}_{\text{tissue}}\;r_{\text{age}} = -0.175$; $\text{VIF}_{B} = 1.412$ vs. $\text{VIF}_{A} = 2.022$)
3. `CST` — **Corticospinal Tract**: Mean of left and right corticospinal tracts ($\text{FA}_{\text{tissue}}\;r_{\text{age}} = -0.111$; $\text{VIF}_{B} = 1.118$)
4. `ATR` — **Anterior Thalamic Radiations**: Operationalized in the JHU ICBM-DTI-81 parcellation as the mean of left and right **anterior corona radiata** ROIs ($\text{FA}_{\text{tissue}}\;r_{\text{age}} = -0.645$; $\text{VIF}_{B} = 2.458$ vs. $\text{VIF}_{A} = 2.885$). _Note_: Retains the conventional ATR label of Kievit et al. (2016), while noting that the volumetric anterior corona radiata ROI encompasses intermingling thalamocortical, corticostriatal, and corticopontine projection fibers.
5. `SLF` — **Superior Longitudinal Fasciculus**: Mean of left and right SLF ($\text{FA}_{\text{tissue}}\;r_{\text{age}} = -0.166$; $\text{VIF}_{B} = 1.394$)
6. `IFOF` — **Inferior Fronto-Occipital Fasciculus**: Mean of left and right IFOF ($\text{FA}_{\text{tissue}}\;r_{\text{age}} = -0.270$; $\text{VIF}_{B} = 1.558$)
7. `ILF` — **Inferior Longitudinal Fasciculus**: Operationalized in the JHU ICBM-DTI-81 parcellation as the mean of left and right **sagittal stratum** ROIs ($\text{FA}_{\text{tissue}}\;r_{\text{age}} = -0.183$; $\text{VIF}_{B} = 1.478$)
8. `UNC` — **Uncinate Fasciculus**: Mean of left and right uncinate fasciculi ($\text{FA}_{\text{tissue}}\;r_{\text{age}} = +0.023$; $\text{VIF}_{B} = 1.273$)
9. `CING` — **Cingulate Cingulum**: Mean of left and right cingulate gyrus cingulum bundles ($\text{FA}_{\text{tissue}}\;r_{\text{age}} = -0.375$; $\text{VIF}_{B} = 1.638$)
10. `CING_Hipp` — **Hippocampal Cingulum**: Mean of left and right parahippocampal cingulum bundles ($\text{FA}_{\text{tissue}}\;r_{\text{age}} = -0.139$; $\text{VIF}_{B} = 1.240$)

### Morphometric & Demographic Covariates

- `age_z` — Standardized chronological age
- `age_sq_z` — Standardized quadratic age $Z(\text{age}_z^2)$
- `sex_z` — Standardized biological sex
- `vent_z` — Standardized log-transformed bilateral lateral ventricular volume normalized by intracranial volume:
  $$\text{vent}_z = Z\!\left(\ln\!\left(\frac{V_{\text{Left Lat Vent}} + V_{\text{Right Lat Vent}}}{\text{eTIV}}\right)\right)$$
