#!/usr/bin/env Rscript
# ==============================================================================
# Script : wa.R (watershed_pipeline.R)
# Purpose: Single-run end-to-end Cam-CAN Watershed Structural Equation Modeling
#          pipeline (QC, feature engineering, JHU tract aggregation, collinearity,
#          4-modality SEM estimation, directionality/mediation tests, OLS
#          covariate residualisation incl. in-memory FreeSurfer lateral
#          ventricle/eTIV adjustment, 3-regime single-pass parallel bootstrap,
#          and outlier sensitivity regimes).
# Usage  : OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript wa.R
#          OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 B_BOOT=5000 Rscript wa.R
# ==============================================================================

# Prevent multi-threaded BLAS contention inside forked mclapply workers
Sys.setenv(
  OMP_NUM_THREADS        = "1",
  OPENBLAS_NUM_THREADS   = "1",
  MKL_NUM_THREADS        = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS    = "1"
)

suppressPackageStartupMessages({
  library(dplyr)
  library(lavaan)
  library(parallel)
})

# ==============================================================================
# 1. CONFIGURATION, PATH RESOLUTION, CONSTANTS & LOGGING ENGINE
# ==============================================================================
input_candidates <- c("output/camcan_watershed_dataset.csv", "camcan_watershed_dataset.csv")
resolved_input   <- input_candidates[file.exists(input_candidates)][1]
if (is.na(resolved_input)) {
  resolved_input <- "output/camcan_watershed_dataset.csv"
}
resolved_out_dir <- dirname(resolved_input)

CFG <- list(
  input_raw = resolved_input,
  out_dir   = resolved_out_dir,
  log_file  = file.path(resolved_out_dir, "watershed_pipeline.log"),
  b_boot    = as.integer(Sys.getenv("B_BOOT", "0")),
  n_cores   = max(1L, detectCores() - 1L),
  seed      = 20261015L
)

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)
if (file.exists(CFG$log_file)) file.remove(CFG$log_file)

RNGkind("L'Ecuyer-CMRG")
set.seed(CFG$seed)

log_msg <- function(level = c("INFO", "WARN", "ERROR", "SUCCESS", "METRIC"), fmt, ...) {
  level <- match.arg(level)
  msg   <- if (length(list(...)) > 0L) sprintf(fmt, ...) else fmt
  line  <- sprintf("[%s] [%-7s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg)
  cat(line, "\n", sep = "")
  cat(line, "\n", file = CFG$log_file, append = TRUE, sep = "")
}

log_section <- function(title) {
  bar   <- strrep("=", 80)
  block <- sprintf("\n%s\n[%s] %s\n%s\n", bar, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), title, bar)
  cat(block, sep = "")
  cat(block, file = CFG$log_file, append = TRUE, sep = "")
}

TRACTS     <- c("FMin", "FMaj", "CST", "ATR", "SLF", "IFOF", "ILF", "UNC", "CING", "CINGHipp")
PS_VARS    <- c("SRTspeed", "SRTcons", "CRTspeed", "CRTcons", "AVspeed", "AVcons")
CAT_VARS   <- paste0("C", 1:4)
MAN_VARS   <- c(CAT_VARS, PS_VARS)
MODALITIES <- c(A = "DTI_FA", B = "FWE_FAtissue", C = "WMTI_AWF", D = "FWE_FW")
MOD_LABELS <- c(
  A = "Model A (Conventional FA)",
  B = "Model B (FWE Tissue FA)",
  C = "Model C (WMTI AWF)",
  D = "Model D (Free-Water FW)"
)
VENT_COLS  <- c("aseg_Left_Lateral_Ventricle_volume", "aseg_Right_Lateral_Ventricle_volume", "eTIV")

tcols <- function(m) paste0(MODALITIES[[m]], "_", TRACTS)
ALL_TRACT_COLS <- unlist(lapply(names(MODALITIES), tcols), use.names = FALSE)

# ==============================================================================
# 2. DATA INGESTION, BEHAVIORAL QC, TRACT AGGREGATION & COHORT PROFILING
# ==============================================================================
log_section("1. DATA INGESTION, QUALITY CONTROL & FEATURE ENGINEERING")
log_msg("INFO", "R: %s | dplyr: %s | lavaan: %s | Cores: %d | B_BOOT: %d",
        R.version.string, as.character(packageVersion("dplyr")),
        as.character(packageVersion("lavaan")), CFG$n_cores, CFG$b_boot)
log_msg("INFO", "Input CSV: %s | Output Dir: %s", CFG$input_raw, CFG$out_dir)

if (!file.exists(CFG$input_raw)) {
  log_msg("ERROR", "Input dataset not found: %s", CFG$input_raw)
  stop(paste("File not found:", CFG$input_raw))
}

raw_df <- read.csv(CFG$input_raw, stringsAsFactors = FALSE, check.names = FALSE)
log_msg("SUCCESS", "Loaded raw dataset: %d rows, %d columns", nrow(raw_df), ncol(raw_df))

# Audit unused extracted columns for physical range anomalies (e.g. wmti_tortuosity overflow)
tort_cols <- grep("^jhu_wmti_tort_", names(raw_df), value = TRUE)
if (length(tort_cols) > 0L) {
  max_tort <- suppressWarnings(max(as.matrix(raw_df[, tort_cols]), na.rm = TRUE))
  if (is.finite(max_tort) && max_tort > 100) {
    log_msg("WARN", "Detected extreme unscaled values in unused jhu_wmti_tort_* columns (max = %.2e). Not used in Models A-D.", max_tort)
  }
}

req_raw_cols <- c(
  "id", "age", "sex",
  "rt_simple_Include_participant", "rt_choice_Include_participant", "mri_behav_include_participant",
  "cattell_SubScore1", "cattell_SubScore2", "cattell_SubScore3", "cattell_SubScore4",
  "rt_simple_RTmean", "rt_simple_RTsd", "rt_choice_RTmean_all", "rt_choice_RTsd_all",
  "mri_behav_mRT", "mri_behav_stdRT",
  VENT_COLS
)
miss_cols <- setdiff(req_raw_cols, names(raw_df))
if (length(miss_cols) > 0L) {
  log_msg("ERROR", "Missing required raw columns: %s", paste(miss_cols, collapse = ", "))
  stop("Missing required columns.")
}

clean_df <- raw_df %>%
  filter(
    !is.na(rt_simple_Include_participant) & rt_simple_Include_participant == 1,
    !is.na(rt_choice_Include_participant) & rt_choice_Include_participant == 1,
    !is.na(mri_behav_include_participant) & mri_behav_include_participant == 1
  ) %>%
  mutate(
    C1        = as.numeric(scale(cattell_SubScore1)),
    C2        = as.numeric(scale(cattell_SubScore2)),
    C3        = as.numeric(scale(cattell_SubScore3)),
    C4        = as.numeric(scale(cattell_SubScore4)),
    SRTspeed  = as.numeric(scale(log(1 / rt_simple_RTmean))),
    SRTcons   = as.numeric(scale(log(1 / rt_simple_RTsd))),
    CRTspeed  = as.numeric(scale(log(1 / rt_choice_RTmean_all))),
    CRTcons   = as.numeric(scale(log(1 / rt_choice_RTsd_all))),
    AV_mean_s = ifelse(mri_behav_mRT > 10, mri_behav_mRT / 1000, mri_behav_mRT),
    AV_sd_s   = ifelse(mri_behav_stdRT > 10, mri_behav_stdRT / 1000, mri_behav_stdRT),
    AVspeed   = as.numeric(scale(log(1 / AV_mean_s))),
    AVcons    = as.numeric(scale(log(1 / AV_sd_s)))
  )

log_msg("INFO", "Behavioral QC: Retained %d / %d participants (%.2f%%)",
        nrow(clean_df), nrow(raw_df), (nrow(clean_df) / nrow(raw_df)) * 100)

get_col <- function(df, col_name) {
  if (!col_name %in% names(df)) {
    log_msg("ERROR", "Missing JHU column: '%s'", col_name)
    stop(paste("Missing column:", col_name))
  }
  df[[col_name]]
}

aggregate_tracts <- function(df, prefix, out_tag) {
  df[[paste0(out_tag, "_FMin")]]     <- as.numeric(scale(get_col(df, paste0(prefix, "Genu of corpus callosum"))))
  df[[paste0(out_tag, "_FMaj")]]     <- as.numeric(scale(get_col(df, paste0(prefix, "Splenium of corpus callosum"))))
  df[[paste0(out_tag, "_CST")]]      <- as.numeric(scale((get_col(df, paste0(prefix, "Corticospinal tract R")) + get_col(df, paste0(prefix, "Corticospinal tract L"))) / 2))
  df[[paste0(out_tag, "_ATR")]]      <- as.numeric(scale((get_col(df, paste0(prefix, "Anterior corona radiata R")) + get_col(df, paste0(prefix, "Anterior corona radiata L"))) / 2))
  df[[paste0(out_tag, "_SLF")]]      <- as.numeric(scale((get_col(df, paste0(prefix, "Superior longitudinal fasciculus R")) + get_col(df, paste0(prefix, "Superior longitudinal fasciculus L"))) / 2))
  df[[paste0(out_tag, "_IFOF")]]     <- as.numeric(scale((get_col(df, paste0(prefix, "Inferior fronto-occipital fasciculus R")) + get_col(df, paste0(prefix, "Inferior fronto-occipital fasciculus L"))) / 2))
  df[[paste0(out_tag, "_ILF")]]      <- as.numeric(scale((get_col(df, paste0(prefix, "Sagittal stratum R")) + get_col(df, paste0(prefix, "Sagittal stratum L"))) / 2))
  df[[paste0(out_tag, "_UNC")]]      <- as.numeric(scale((get_col(df, paste0(prefix, "Uncinate fasciculus R")) + get_col(df, paste0(prefix, "Uncinate fasciculus L"))) / 2))
  df[[paste0(out_tag, "_CING")]]     <- as.numeric(scale((get_col(df, paste0(prefix, "Cingulate cingulum R")) + get_col(df, paste0(prefix, "Cingulate cingulum L"))) / 2))
  df[[paste0(out_tag, "_CINGHipp")]] <- as.numeric(scale((get_col(df, paste0(prefix, "Hippocampal cingulum R")) + get_col(df, paste0(prefix, "Hippocampal cingulum L"))) / 2))
  df
}

clean_df <- aggregate_tracts(clean_df, "jhu_dti_fa_",   "DTI_FA")
clean_df <- aggregate_tracts(clean_df, "jhu_fwe_fa_",   "FWE_FAtissue")
clean_df <- aggregate_tracts(clean_df, "jhu_wmti_awf_", "WMTI_AWF")
clean_df <- aggregate_tracts(clean_df, "jhu_fwe_frac_", "FWE_FW")

clean_df$outlier_flag_behav   <- apply(clean_df[, MAN_VARS], 1, function(r) any(abs(r) > 4, na.rm = TRUE))
clean_df$outlier_flag_tract   <- apply(clean_df[, ALL_TRACT_COLS], 1, function(r) any(abs(r) > 4, na.rm = TRUE))
clean_df$outlier_flag         <- clean_df$outlier_flag_behav | clean_df$outlier_flag_tract
clean_df$complete_behav_flag  <- apply(clean_df[, MAN_VARS], 1, function(r) !any(is.na(r)))
clean_df$complete_dwi_flag    <- apply(clean_df[, ALL_TRACT_COLS], 1, function(r) !any(is.na(r)))
clean_df$analysis_cohort_flag <- (!clean_df$outlier_flag) & clean_df$complete_behav_flag & clean_df$complete_dwi_flag

id_col  <- "id"
sex_col <- "sex"

analysis_df <- clean_df %>% filter(analysis_cohort_flag == TRUE)
log_msg("SUCCESS", "Harmonized Analysis Cohort: N = %d (Outliers: %d | Incomplete DWI: %d | Incomplete Behav: %d)",
        nrow(analysis_df), sum(clean_df$outlier_flag), sum(!clean_df$complete_dwi_flag), sum(!clean_df$complete_behav_flag))

write.csv(clean_df, file.path(CFG$out_dir, "watershed_sem_ready.csv"), row.names = FALSE)

all_manifest <- c(MAN_VARS, ALL_TRACT_COLS)
preproc_summary <- data.frame(
  Feature       = all_manifest,
  Mean          = round(sapply(all_manifest, function(v) mean(clean_df[[v]], na.rm = TRUE)), 4),
  SD            = round(sapply(all_manifest, function(v) sd(clean_df[[v]], na.rm = TRUE)), 4),
  Min           = round(sapply(all_manifest, function(v) min(clean_df[[v]], na.rm = TRUE)), 4),
  Max           = round(sapply(all_manifest, function(v) max(clean_df[[v]], na.rm = TRUE)), 4),
  Missing_Count = sapply(all_manifest, function(v) sum(is.na(clean_df[[v]]))),
  row.names     = NULL
)
write.csv(preproc_summary, file.path(CFG$out_dir, "preprocessing_summary.csv"), row.names = FALSE)

# Demographics, Age Correlations & Attrition Profile
q_age <- quantile(analysis_df$age, c(0.25, 0.50, 0.75), na.rm = TRUE)
demo_df <- data.frame(
  n          = nrow(analysis_df),
  age_mean   = round(mean(analysis_df$age, na.rm = TRUE), 2),
  age_sd     = round(sd(analysis_df$age, na.rm = TRUE), 2),
  age_min    = round(min(analysis_df$age, na.rm = TRUE), 1),
  age_max    = round(max(analysis_df$age, na.rm = TRUE), 1),
  age_median = round(q_age[2], 2),
  age_q25    = round(q_age[1], 2),
  age_q75    = round(q_age[3], 2),
  row.names  = NULL
)
write.csv(demo_df, file.path(CFG$out_dir, "cohort_demographics.csv"), row.names = FALSE)
write.csv(as.data.frame(table(analysis_df[[sex_col]])), file.path(CFG$out_dir, "cohort_sex_counts.csv"), row.names = FALSE)

agecor_df <- data.frame(
  variable       = all_manifest,
  r_age_pearson  = round(sapply(all_manifest, function(v) cor(analysis_df$age, analysis_df[[v]], use = "complete.obs", method = "pearson")), 3),
  r_age_spearman = round(sapply(all_manifest, function(v) cor(analysis_df$age, analysis_df[[v]], use = "complete.obs", method = "spearman")), 3),
  row.names      = NULL
)
write.csv(agecor_df, file.path(CFG$out_dir, "age_correlations.csv"), row.names = FALSE)

ex_df <- clean_df[!clean_df$analysis_cohort_flag, ]
if (nrow(ex_df) > 0L) {
  prof_df <- data.frame(
    id                = ex_df[[id_col]],
    age               = ex_df$age,
    sex               = ex_df[[sex_col]],
    behav_outlier     = ex_df$outlier_flag_behav,
    tract_outlier     = ex_df$outlier_flag_tract,
    dwi_incomplete    = !ex_df$complete_dwi_flag,
    flagged_variables = apply(ex_df[, all_manifest, drop = FALSE], 1, function(r) paste(names(r)[which(abs(r) > 4)], collapse = ";")),
    max_abs_z         = round(apply(ex_df[, all_manifest, drop = FALSE], 1, function(r) {
      vals <- abs(r[!is.na(r)])
      if (length(vals) == 0L) NA_real_ else max(vals)
    }), 2),
    row.names         = NULL
  )
} else {
  prof_df <- data.frame(
    id = character(), age = numeric(), sex = character(),
    behav_outlier = logical(), tract_outlier = logical(), dwi_incomplete = logical(),
    flagged_variables = character(), max_abs_z = numeric()
  )
}
write.csv(prof_df, file.path(CFG$out_dir, "excluded_participants_profile.csv"), row.names = FALSE)

safe_wilcox <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) < 1L || length(y) < 1L) return(NA_real_)
  tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
}

cmp_age_df <- data.frame(
  group    = c("retained", "excluded_all", "excluded_tract_outlier", "excluded_dwi_incomplete"),
  n        = c(nrow(analysis_df), nrow(ex_df), sum(ex_df$outlier_flag_tract), sum(!ex_df$complete_dwi_flag)),
  mean_age = round(c(
    mean(analysis_df$age, na.rm = TRUE),
    if (nrow(ex_df) > 0L) mean(ex_df$age, na.rm = TRUE) else NA_real_,
    if (sum(ex_df$outlier_flag_tract) > 0L) mean(ex_df$age[ex_df$outlier_flag_tract], na.rm = TRUE) else NA_real_,
    if (sum(!ex_df$complete_dwi_flag) > 0L) mean(ex_df$age[!ex_df$complete_dwi_flag], na.rm = TRUE) else NA_real_
  ), 2),
  wilcox_p_vs_retained = signif(c(
    NA_real_,
    safe_wilcox(ex_df$age, analysis_df$age),
    safe_wilcox(ex_df$age[ex_df$outlier_flag_tract], analysis_df$age),
    safe_wilcox(ex_df$age[!ex_df$complete_dwi_flag], analysis_df$age)
  ), 3),
  row.names = NULL
)
write.csv(cmp_age_df, file.path(CFG$out_dir, "excluded_vs_retained_age.csv"), row.names = FALSE)

# ==============================================================================
# 3. COLLINEARITY DIAGNOSTICS & SUMMARY DEPOSIT MATRICES
# ==============================================================================
log_section("2. COLLINEARITY DIAGNOSTICS & SUMMARY COVARIANCE DEPOSITS")

vif_list <- list()
for (m in names(MODALITIES)) {
  v_mod <- c(MAN_VARS, tcols(m))
  write.csv(cov(analysis_df[, v_mod], use = "complete.obs"),
            file.path(CFG$out_dir, sprintf("covariance_matrix_%s.csv", m)))
  write.csv(data.frame(variable = v_mod, mean = colMeans(analysis_df[, v_mod], na.rm = TRUE), row.names = NULL),
            file.path(CFG$out_dir, sprintf("means_%s.csv", m)), row.names = FALSE)

  R_mat <- cor(analysis_df[, tcols(m)], use = "complete.obs")
  colnames(R_mat) <- rownames(R_mat) <- TRACTS
  write.csv(round(R_mat, 3), file.path(CFG$out_dir, sprintf("tract_correlations_%s.csv", m)))
  vif_list[[m]] <- data.frame(model = paste0("Model_", m), modality = MODALITIES[[m]],
                              tract = TRACTS, VIF = round(diag(solve(R_mat)), 3), row.names = NULL)
}
write.csv(do.call(rbind, vif_list), file.path(CFG$out_dir, "tract_vif_all_models.csv"), row.names = FALSE)
log_msg("SUCCESS", "Exported tract correlation matrices, VIFs, and deposit covariance/mean tables.")

# ==============================================================================
# 4. UNIFIED SEM SYNTAX GENERATOR & EXTRACTION HELPERS
# ==============================================================================
build_ws_syntax <- function(prefix, kind = c("forward", "direct", "indirect", "reversed", "inverted", "wmonly")) {
  kind    <- match.arg(kind)
  T_vars  <- paste0(prefix, "_", TRACTS)
  meas    <- "FI =~ C1 + C2 + C3 + C4"
  ps_cov  <- vapply(1:5, function(i) sprintf("%s ~~ %s", PS_VARS[i], paste(PS_VARS[(i + 1):6], collapse = " + ")), "")
  t_cov   <- vapply(1:9, function(i) sprintf("%s ~~ %s", T_vars[i], paste(T_vars[(i + 1):10], collapse = " + ")), "")
  fi_ps   <- paste0("FI ~ ", paste(sprintf("f_%s*%s", PS_VARS, PS_VARS), collapse = " + "))
  ps_reg  <- vapply(PS_VARS, function(o) {
    sprintf("%s ~ %s", o, paste(sprintf("%s*%s", paste0("p_", o, "_", TRACTS), T_vars), collapse = " + "))
  }, "")
  rev_reg <- vapply(T_vars, function(t) sprintf("%s ~ %s", t, paste(PS_VARS, collapse = " + ")), "")
  fi_wm   <- paste0("FI ~ ", paste(sprintf("d_%s*%s", TRACTS, T_vars), collapse = " + "))
  ind_def <- vapply(TRACTS, function(t) {
    sprintf("ind_%s := %s", t, paste(sprintf("p_%s_%s*f_%s", PS_VARS, t, PS_VARS), collapse = " + "))
  }, "")

  parts <- switch(
    kind,
    forward  = c(meas, fi_ps, ps_reg, ps_cov, t_cov),
    direct   = c(meas, fi_ps, fi_wm, ps_reg, ps_cov, t_cov),
    indirect = c(meas, fi_ps, ps_reg, ps_cov, t_cov, ind_def),
    reversed = c(meas, rev_reg, paste("FI ~", paste(PS_VARS, collapse = " + ")), paste0("FI ~~ 0*", T_vars), ps_cov, t_cov),
    inverted = c(meas, rev_reg, paste("FI ~", paste(T_vars, collapse = " + ")), ps_cov, t_cov),
    wmonly   = c(meas, paste("FI ~", paste(T_vars, collapse = " + ")), t_cov)
  )
  paste(parts, collapse = "\n")
}

WS_FORWARD_SYNTAX <- setNames(
  lapply(names(MODALITIES), function(m) build_ws_syntax(MODALITIES[[m]], "forward")),
  names(MODALITIES)
)

fit_ws <- function(m, kind = "forward", data = analysis_df) {
  suppressWarnings(
    sem(build_ws_syntax(MODALITIES[[m]], kind), data = data, estimator = "MLR", missing = "fiml", fixed.x = FALSE)
  )
}

# Fast point-estimate fitter for bootstrap replicates (identical ML objective, no SE/Jacobian overhead)
fit_ws_boot <- function(m, data) {
  has_na <- anyNA(data[, c(MAN_VARS, tcols(m))])
  suppressWarnings(
    sem(
      WS_FORWARD_SYNTAX[[m]],
      data          = data,
      estimator     = "ML",
      likelihood    = "normal",
      meanstructure = TRUE,
      missing       = if (has_na) "fiml" else "listwise",
      fixed.x       = FALSE,
      se            = "none",
      test          = "none"
    )
  )
}

fit_row <- function(fit, name) {
  nm <- c("chisq.scaled", "df", "pvalue.scaled", "cfi.robust", "tli.robust",
          "rmsea.robust", "rmsea.ci.lower.robust", "rmsea.ci.upper.robust", "srmr", "aic", "bic",
          "scaling.factor.h0", "scaling.factor.h1")
  fm <- tryCatch(fitMeasures(fit, nm), error = function(e) setNames(rep(NA_real_, length(nm)), nm))
  vals <- as.list(setNames(round(as.numeric(fm), 4), c(
    "chisq.scaled", "df", "pvalue.scaled", "cfi.robust", "tli.robust",
    "rmsea.robust", "rmsea.ci.lower.robust", "rmsea.ci.upper.robust", "srmr", "aic", "bic", "c_H0", "c_H1"
  )))
  data.frame(
    model     = name,
    n         = lavInspect(fit, "nobs"),
    converged = lavInspect(fit, "converged"),
    as.data.frame(vals),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

paths_tbl <- function(fit, m, label) {
  s <- standardizedSolution(fit)
  s <- s[s$op == "~" & s$lhs %in% PS_VARS & s$rhs %in% tcols(m), ]
  s$tract  <- sub(paste0("^", MODALITIES[[m]], "_"), "", s$rhs)
  s$model  <- label
  s$p_holm <- p.adjust(s$pvalue, "holm")
  s$p_bh   <- p.adjust(s$pvalue, "BH")
  s[, c("model", "lhs", "tract", "est.std", "se", "z", "pvalue", "ci.lower", "ci.upper", "p_holm", "p_bh")]
}

# Scale-aware Wald test: tests equality of standardized betas (b_ATR * SD_ATR == b_FMin * SD_FMin)
extract_wald_tbl <- function(fits_list, model_keys = c("A", "B"), prefix_tag = "Model", holm_by_model = TRUE) {
  w_df <- do.call(rbind, lapply(model_keys, function(m) {
    cov_ov  <- lavInspect(fits_list[[m]], "cov.ov")
    atr_var <- paste0(MODALITIES[[m]], "_ATR")
    fmn_var <- paste0(MODALITIES[[m]], "_FMin")
    sd_atr  <- sqrt(as.numeric(cov_ov[atr_var, atr_var]))
    sd_fmn  <- sqrt(as.numeric(cov_ov[fmn_var, fmn_var]))

    do.call(rbind, lapply(PS_VARS, function(o) {
      con_str <- sprintf("%.12f * p_%s_ATR == %.12f * p_%s_FMin", sd_atr, o, sd_fmn, o)
      w <- lavTestWald(fits_list[[m]], con_str)
      data.frame(
        model      = paste0(prefix_tag, "_", m),
        outcome    = o,
        wald_chisq = round(as.numeric(w$stat), 3),
        df         = w$df,
        p          = signif(as.numeric(w$p.value), 4),
        row.names  = NULL
      )
    }))
  }))
  w_df$p_holm <- if (holm_by_model) {
    ave(w_df$p, w_df$model, FUN = function(p) p.adjust(p, "holm"))
  } else {
    p.adjust(w_df$p, "holm")
  }
  w_df
}

# ==============================================================================
# 5. PRIMARY WATERSHED MODELS (A-D), STANDALONE CFA & WALD TESTS
# ==============================================================================
log_section("3. ESTIMATING PRIMARY WATERSHED MODELS (A-D) & MEASUREMENT CFA")

cfa_fit <- suppressWarnings(cfa("FI =~ C1 + C2 + C3 + C4", data = analysis_df, estimator = "MLR", missing = "fiml"))
write.csv(fit_row(cfa_fit, "CFA_FI_alone"), file.path(CFG$out_dir, "cfa_fi_fit.csv"), row.names = FALSE)

base_fits <- setNames(lapply(names(MODALITIES), function(m) {
  f <- fit_ws(m, "forward", analysis_df)
  log_msg("SUCCESS", "%s | ChiSq = %.3f (df = %.0f), CFI = %.4f, RMSEA = %.4f, SRMR = %.4f, AIC = %.2f",
          MOD_LABELS[[m]], fitMeasures(f, "chisq.scaled"), fitMeasures(f, "df"),
          fitMeasures(f, "cfi.robust"), fitMeasures(f, "rmsea.robust"), fitMeasures(f, "srmr"), AIC(f))
  f
}), names(MODALITIES))

model_comp_df <- do.call(rbind, Map(fit_row, base_fits, unname(MOD_LABELS[names(base_fits)]))) %>%
  arrange(rmsea.robust, desc(cfi.robust)) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, everything())
write.csv(model_comp_df, file.path(CFG$out_dir, "watershed_model_comparison.csv"), row.names = FALSE)

write.csv(
  do.call(rbind, Map(function(f, m) paths_tbl(f, m, paste0("Model_", m)), base_fits, names(base_fits))),
  file.path(CFG$out_dir, "tract_to_speed_paths_unadjusted.csv"),
  row.names = FALSE
)

std_all <- do.call(rbind, lapply(names(base_fits), function(m) {
  s <- standardizedSolution(base_fits[[m]])
  s$model <- paste0("Model_", m)
  s
}))
unstd_all <- do.call(rbind, lapply(names(base_fits), function(m) {
  u <- parameterEstimates(base_fits[[m]])
  u$model <- paste0("Model_", m)
  u
}))
write.csv(std_all,   file.path(CFG$out_dir, "standardized_solutions_all_models.csv"), row.names = FALSE)
write.csv(unstd_all, file.path(CFG$out_dir, "unstandardized_solutions_all_models.csv"), row.names = FALSE)

r2_unadj <- do.call(rbind, lapply(names(base_fits), function(m) {
  r <- lavInspect(base_fits[[m]], "rsquare")[c(PS_VARS, "FI")]
  data.frame(model = paste0("Model_", m), var = names(r), r2 = round(r, 4), row.names = NULL)
}))
write.csv(r2_unadj, file.path(CFG$out_dir, "rsquared_unadjusted.csv"), row.names = FALSE)

wald_unadj <- extract_wald_tbl(base_fits, c("A", "B"), prefix_tag = "Model", holm_by_model = TRUE)
write.csv(wald_unadj, file.path(CFG$out_dir, "wald_ATR_vs_FMin_unadjusted.csv"), row.names = FALSE)

# ==============================================================================
# 6. DIRECTIONALITY, HIERARCHY INVERSION & MEDIATION ANALYSES
# ==============================================================================
log_section("4. DIRECTIONALITY, HIERARCHY INVERSION & MEDIATION TESTS")

dir_rows  <- list()
lrt_rows  <- list()
dir_paths <- list()
ind_rows  <- list()
r2fi_rows <- list()

for (m in c("A", "B")) {
  F0 <- base_fits[[m]]
  Fd <- fit_ws(m, "direct", analysis_df)
  E  <- fit_ws(m, "reversed", analysis_df)
  I  <- fit_ws(m, "inverted", analysis_df)
  Fi <- fit_ws(m, "indirect", analysis_df)
  Fw <- fit_ws(m, "wmonly", analysis_df)

  d_sub <- rbind(
    fit_row(F0, paste0("Model_", m, "_forward")),
    fit_row(Fd, paste0("Model_", m, "_forward_plus_direct_WM_to_FI")),
    fit_row(E,  paste0("Model_", m, "_reversed_WM_PS_arrows")),
    fit_row(I,  paste0("Model_", m, "_inverted_hierarchy"))
  )
  d_sub$delta_aic_vs_forward <- round(d_sub$aic - AIC(F0), 2)
  dir_rows[[m]] <- d_sub

  a_lrt <- lavTestLRT(F0, Fd)
  lrt_rows[[m]] <- data.frame(
    model                         = paste0("Model_", m),
    test                          = "forward vs forward+10 direct WM->FI (scaled LRT)",
    chisq_diff                    = round(a_lrt[2, "Chisq diff"], 3),
    df_diff                       = a_lrt[2, "Df diff"],
    p                             = signif(a_lrt[2, "Pr(>Chisq)"], 4),
    logLik_forward_minus_reversed = signif(as.numeric(logLik(F0) - logLik(E)), 4),
    df_forward                    = fitMeasures(F0, "df"),
    df_reversed                   = fitMeasures(E, "df"),
    delta_aic_inverted_vs_forward = round(AIC(I) - AIC(F0), 2),
    row.names                     = NULL
  )

  s_dir <- standardizedSolution(Fd)
  s_dir <- s_dir[s_dir$op == "~" & s_dir$lhs == "FI" & s_dir$rhs %in% tcols(m), ]
  s_dir$model <- paste0("Model_", m)
  dir_paths[[m]] <- s_dir

  pe_unstd <- parameterEstimates(Fi, standardized = TRUE)
  pe_unstd <- pe_unstd[pe_unstd$op == ":=", ]
  pe_std   <- standardizedSolution(Fi)
  pe_std   <- pe_std[pe_std$op == ":=", ]

  ind_rows[[m]] <- data.frame(
    lhs            = pe_std$lhs,
    est            = pe_unstd$est,
    se             = pe_std$se,
    z              = pe_std$z,
    pvalue         = pe_std$pvalue,
    ci.lower       = pe_std$ci.lower,
    ci.upper       = pe_std$ci.upper,
    std.all        = pe_std$est.std,
    se_unstd       = pe_unstd$se,
    ci.lower_unstd = pe_unstd$ci.lower,
    ci.upper_unstd = pe_unstd$ci.upper,
    model          = paste0("Model_", m),
    p_holm         = p.adjust(pe_std$pvalue, "holm"),
    row.names      = NULL,
    stringsAsFactors = FALSE
  )

  r2fi_rows[[m]] <- data.frame(
    model                     = paste0("Model_", m),
    R2_FI_from_PS_via_forward = round(lavInspect(F0, "rsquare")["FI"], 4),
    R2_FI_from_tracts_only    = round(lavInspect(Fw, "rsquare")["FI"], 4),
    R2_FI_forward_plus_direct = round(lavInspect(Fd, "rsquare")["FI"], 4),
    row.names                 = NULL
  )
}

write.csv(do.call(rbind, dir_rows),  file.path(CFG$out_dir, "directionality_model_comparison.csv"), row.names = FALSE)
write.csv(do.call(rbind, lrt_rows),  file.path(CFG$out_dir, "mediation_and_equivalence_tests.csv"), row.names = FALSE)
write.csv(do.call(rbind, dir_paths), file.path(CFG$out_dir, "direct_wm_to_fi_paths.csv"), row.names = FALSE)
write.csv(do.call(rbind, ind_rows),  file.path(CFG$out_dir, "indirect_effects_delta_method.csv"), row.names = FALSE)
write.csv(do.call(rbind, r2fi_rows), file.path(CFG$out_dir, "fi_rsquared_decomposition.csv"), row.names = FALSE)
log_msg("SUCCESS", "Completed hierarchy inversion, LRT mediation, and indirect effect decompositions.")

# ==============================================================================
# 7. COVARIATE-RESIDUALISED WATERSHED MODELS (df = 60)
# ==============================================================================
log_section("5. COVARIATE-RESIDUALISED WATERSHED MODELS (AGE, AGE^2, SEX, VENTRICLE/eTIV)")

add_covars <- function(d) {
  d$age_z  <- as.numeric(scale(d$age))
  d$age2_z <- as.numeric(scale(d$age_z^2))
  cv <- c("age_z", "age2_z")
  if (!is.na(sex_col) && sex_col %in% names(d)) {
    d$sex_z <- as.numeric(scale(as.numeric(factor(d[[sex_col]]))))
    cv <- c(cv, "sex_z")
  }
  if ("motion_rms" %in% names(d) && !all(is.na(d$motion_rms))) {
    d$motion_z <- as.numeric(scale(d$motion_rms))
    cv <- c(cv, "motion_z")
  }
  if (all(VENT_COLS %in% names(d))) {
    vent_icv <- (d$aseg_Left_Lateral_Ventricle_volume + d$aseg_Right_Lateral_Ventricle_volume) / d$eTIV
    d$vent_z <- as.numeric(scale(log(vent_icv)))
  }
  list(d = d, cv = cv)
}

# Vectorized OLS residualisation with z-standardisation of residuals
residualise <- function(d, cv, vars) {
  sub_vars <- d[, vars, drop = FALSE]
  sub_cv   <- d[, cv, drop = FALSE]
  if (!anyNA(sub_vars) && !anyNA(sub_cv)) {
    X   <- model.matrix(~ ., data = sub_cv)
    res <- qr.resid(qr(X), as.matrix(sub_vars))
    d[, vars] <- lapply(seq_along(vars), function(j) as.numeric(scale(res[, j])))
  } else {
    f <- as.formula(paste("y ~", paste(cv, collapse = " + ")))
    for (v in vars) {
      tmp    <- data.frame(y = d[[v]], sub_cv)
      d[[v]] <- as.numeric(scale(resid(lm(f, tmp, na.action = na.exclude))))
    }
  }
  d
}

run_adjusted_suite <- function(d0, tag, extra_cv = character()) {
  ac <- add_covars(d0)
  cv <- c(ac$cv, extra_cv)
  dr <- residualise(ac$d, cv, c(MAN_VARS, ALL_TRACT_COLS))
  log_msg("INFO", "[%s] Fitting residualised SEM | Covariates: %s | N = %d", tag, paste(cv, collapse = ", "), nrow(dr))
  fits <- setNames(lapply(names(MODALITIES), function(m) fit_ws(m, "forward", dr)), names(MODALITIES))

  r2_df <- do.call(rbind, lapply(names(fits), function(m) {
    r <- lavInspect(fits[[m]], "rsquare")[c(PS_VARS, "FI")]
    data.frame(model = paste0(tag, "_Model_", m), var = names(r), r2 = round(r, 4), row.names = NULL)
  }))

  list(
    fits  = fits,
    fit   = do.call(rbind, Map(fit_row, fits, paste0(tag, "_Model_", names(fits)))),
    paths = do.call(rbind, Map(function(f, m) paths_tbl(f, m, paste0(tag, "_Model_", m)), fits, names(fits))),
    wald  = extract_wald_tbl(fits, c("A", "B"), prefix_tag = paste0(tag, "_Model"), holm_by_model = TRUE),
    r2    = r2_df
  )
}

adj_res <- run_adjusted_suite(analysis_df, "adj_age_sex")
adj_fits_all  <- adj_res$fit
adj_paths_all <- adj_res$paths
adj_wald_all  <- adj_res$wald
adj_r2_all    <- adj_res$r2

vent_res <- NULL
if (all(VENT_COLS %in% names(analysis_df))) {
  d_vent <- analysis_df %>%
    filter(!is.na(aseg_Left_Lateral_Ventricle_volume),
           !is.na(aseg_Right_Lateral_Ventricle_volume),
           !is.na(eTIV), eTIV > 0)
  vent_res      <- run_adjusted_suite(d_vent, "adj_age_sex_vent", "vent_z")
  adj_fits_all  <- rbind(adj_fits_all,  vent_res$fit)
  adj_paths_all <- rbind(adj_paths_all, vent_res$paths)
  adj_wald_all  <- rbind(adj_wald_all,  vent_res$wald)
  adj_r2_all    <- rbind(adj_r2_all,    vent_res$r2)
  log_msg("SUCCESS", "Fitted ventricle-adjusted models from in-memory aseg_volume + eTIV (N = %d).", nrow(d_vent))
} else {
  log_msg("WARN", "Required aseg lateral ventricle columns not found; skipping ventricle model.")
}

write.csv(adj_fits_all,  file.path(CFG$out_dir, "watershed_adjusted_fit.csv"), row.names = FALSE)
write.csv(adj_paths_all, file.path(CFG$out_dir, "tract_to_speed_paths_adjusted.csv"), row.names = FALSE)
write.csv(adj_wald_all,  file.path(CFG$out_dir, "wald_ATR_vs_FMin_adjusted.csv"), row.names = FALSE)
write.csv(adj_r2_all,    file.path(CFG$out_dir, "rsquared_adjusted.csv"), row.names = FALSE)
log_msg("SUCCESS", "Exported covariate-residualised model fits, paths, Wald tests, and R2.")

# ==============================================================================
# 8. SINGLE-PASS PARALLEL BOOTSTRAP (UNADJUSTED, AGE/SEX & VENTRICLE-ADJUSTED)
# ==============================================================================
log_section("6. SINGLE-PASS PARALLEL BOOTSTRAP OF ATR vs. FMin CONTRASTS")

get_focal_betas <- function(fit, m, tag) {
  s <- standardizedSolution(fit, se = FALSE, zstat = FALSE, pvalue = FALSE, ci = FALSE)
  s <- s[s$op == "~" & s$lhs %in% PS_VARS & s$rhs %in% paste0(MODALITIES[[m]], "_", c("FMin", "ATR")), ]
  setNames(s$est.std, paste(tag, m, s$lhs, sub(paste0("^", MODALITIES[[m]], "_"), "", s$rhs), sep = "|"))
}

focal_boot_vars <- c(MAN_VARS, tcols("A"), tcols("B"))
has_vent_boot   <- !is.null(vent_res)

boot_worker <- function(i) {
  d_unadj <- analysis_df[sample.int(nrow(analysis_df), replace = TRUE), ]
  ac      <- add_covars(d_unadj)
  d_adj   <- residualise(ac$d, ac$cv, focal_boot_vars)
  d_vent  <- if (has_vent_boot && "vent_z" %in% names(ac$d)) {
    residualise(ac$d, c(ac$cv, "vent_z"), focal_boot_vars)
  } else {
    NULL
  }

  tryCatch({
    fa_u <- fit_ws_boot("A", d_unadj)
    fb_u <- fit_ws_boot("B", d_unadj)
    fa_a <- fit_ws_boot("A", d_adj)
    fb_a <- fit_ws_boot("B", d_adj)
    if (!lavInspect(fa_u, "converged") || !lavInspect(fb_u, "converged") ||
        !lavInspect(fa_a, "converged") || !lavInspect(fb_a, "converged")) {
      return(NULL)
    }
    out <- c(
      get_focal_betas(fa_u, "A", "unadj"),
      get_focal_betas(fb_u, "B", "unadj"),
      get_focal_betas(fa_a, "A", "adj"),
      get_focal_betas(fb_a, "B", "adj")
    )
    if (!is.null(d_vent)) {
      fa_v <- fit_ws_boot("A", d_vent)
      fb_v <- fit_ws_boot("B", d_vent)
      if (!lavInspect(fa_v, "converged") || !lavInspect(fb_v, "converged")) {
        return(NULL)
      }
      out <- c(
        out,
        get_focal_betas(fa_v, "A", "adj_vent"),
        get_focal_betas(fb_v, "B", "adj_vent")
      )
    }
    out
  }, error = function(e) NULL)
}

log_msg("INFO", "Executing single-pass parallel bootstrap (B = %d, Cores = %d)...", CFG$b_boot, CFG$n_cores)
boot_list <- mclapply(seq_len(CFG$b_boot), boot_worker, mc.cores = CFG$n_cores, mc.set.seed = TRUE)
boot_mat  <- do.call(rbind, Filter(Negate(is.null), boot_list))
log_msg("SUCCESS", "Bootstrap converged across all specifications: %d / %d replicates", nrow(boot_mat), CFG$b_boot)

summarize_boot <- function(tag, ref_fits, spec_label) {
  ref_std <- lapply(ref_fits[c("A", "B")], function(f) {
    standardizedSolution(f, se = FALSE, zstat = FALSE, pvalue = FALSE, ci = FALSE)
  })
  obs_beta <- function(m, o, t) {
    s <- ref_std[[m]]
    s$est.std[s$op == "~" & s$lhs == o & s$rhs == paste0(MODALITIES[[m]], "_", t)]
  }
  rows <- list()
  for (o in PS_VARS) {
    g <- function(m, t) boot_mat[, paste(tag, m, o, t, sep = "|")]
    defs <- list(
      ATR_minus_FMin_in_A       = list(g("A", "ATR") - g("A", "FMin"), obs_beta("A", o, "ATR") - obs_beta("A", o, "FMin")),
      ATR_minus_FMin_in_B       = list(g("B", "ATR") - g("B", "FMin"), obs_beta("B", o, "ATR") - obs_beta("B", o, "FMin")),
      change_in_contrast_B_vs_A = list((g("B", "ATR") - g("B", "FMin")) - (g("A", "ATR") - g("A", "FMin")),
                                       (obs_beta("B", o, "ATR") - obs_beta("B", o, "FMin")) - (obs_beta("A", o, "ATR") - obs_beta("A", o, "FMin"))),
      FMin_B_minus_A            = list(g("B", "FMin") - g("A", "FMin"), obs_beta("B", o, "FMin") - obs_beta("A", o, "FMin")),
      ATR_B_minus_A             = list(g("B", "ATR")  - g("A", "ATR"),  obs_beta("B", o, "ATR")  - obs_beta("A", o, "ATR"))
    )
    for (k in names(defs)) {
      rows[[length(rows) + 1L]] <- data.frame(
        specification = spec_label,
        outcome       = o,
        contrast      = k,
        estimate      = round(defs[[k]][[2]], 4),
        ci_lo         = round(unname(quantile(defs[[k]][[1]], 0.025, na.rm = TRUE)), 4),
        ci_hi         = round(unname(quantile(defs[[k]][[1]], 0.975, na.rm = TRUE)), 4),
        n_boot        = nrow(boot_mat),
        row.names     = NULL
      )
    }
  }
  do.call(rbind, rows)
}

boot_summary_df <- rbind(
  summarize_boot("unadj", base_fits,    "unadjusted"),
  summarize_boot("adj",   adj_res$fits, "covariate_residualised_age_sex")
)
if (has_vent_boot) {
  boot_summary_df <- rbind(
    boot_summary_df,
    summarize_boot("adj_vent", vent_res$fits, "covariate_residualised_age_sex_vent")
  )
}
write.csv(boot_summary_df, file.path(CFG$out_dir, "bootstrap_ATR_FMin_contrasts.csv"), row.names = FALSE)
log_msg("SUCCESS", "Exported unified bootstrap contrasts table.")

# ==============================================================================
# 9. UNIFIED OUTLIER & SAMPLE-PARITY SENSITIVITY MODULE
# ==============================================================================
log_section("7. UNIFIED OUTLIER & ATTRITION SENSITIVITY ANALYSIS")

sens_fits  <- list()
sens_paths <- list()
common_noexcl_df <- clean_df[complete.cases(clean_df[, all_manifest]), ]

for (sc in c("no_outlier_exclusion", "modality_specific_4SD", "common_noexcl")) {
  for (m in names(MODALITIES)) {
    vars_m <- c(MAN_VARS, tcols(m))
    d_sub <- switch(
      sc,
      no_outlier_exclusion  = clean_df[complete.cases(clean_df[, vars_m]), ],
      modality_specific_4SD = {
        d_cc <- clean_df[complete.cases(clean_df[, vars_m]), ]
        d_cc[apply(d_cc[, vars_m], 1, function(r) all(abs(r) <= 4)), ]
      },
      common_noexcl         = common_noexcl_df
    )
    f_sub <- fit_ws(m, "forward", d_sub)
    lab   <- paste0(sc, "_Model_", m)
    sens_fits[[lab]]  <- fit_row(f_sub, lab)
    sens_paths[[lab]] <- paths_tbl(f_sub, m, lab)
  }
}

write.csv(do.call(rbind, sens_fits),  file.path(CFG$out_dir, "sensitivity_outlier_regimes_fit.csv"), row.names = FALSE)
write.csv(do.call(rbind, sens_paths), file.path(CFG$out_dir, "sensitivity_outlier_regimes_paths.csv"), row.names = FALSE)
log_msg("SUCCESS", "Exported unified sensitivity tables across all 3 outlier regimes.")

# ==============================================================================
# 10. SESSION INFO & COMPLETION
# ==============================================================================
writeLines(capture.output(sessionInfo()), file.path(CFG$out_dir, "sessionInfo.txt"))
log_section("WATERSHED SEM SINGLE-RUN PIPELINE COMPLETED SUCCESSFULLY")
