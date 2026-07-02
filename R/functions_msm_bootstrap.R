################################################################################
# functions_msm_bootstrap.R
#
# Shared functions for the marginal structural models (MSM), the nonparametric
# bootstrap (stratified cluster resampling of individuals), and extraction of
# 5-year risks / risk differences with bootstrap confidence intervals.
#
# The MSM differs by trial, so each trial has its own outcome-model function:
#   - Trial 1: analyze_survival_outcomes_trial1()  (categorical arm)
#   - Trial 2: analyze_survival_outcomes_trial2()  (dose-response for the nine
#              interruption windows + a categorical model for the controls)
# Both are covariate-adjusted MSMs; strategy-specific cumulative incidences are
# obtained by g-computation (predict for every individual under each strategy,
# then average over the baseline covariate distribution).
#
# calculate_ip_weights() is defined in the trial's CCW script (01_/02_).
# bootstrap_survival_analysis() takes the trial's MSM function via `analyze_fun`.
################################################################################

# Baseline covariates included as main effects in every outcome model
BASE_COVARS <- paste(
  "age_cl", "previous_pregnancy", "breast_surgery_2cl", "ct_setting_5cl", "rt",
  "ht_type", "pnuicc_2cl_modified", "her2_status", "nb_comor_ql", "Cardiovascular",
  "Endocrine_and_metabolism", "Psychiatric_disorders", "Other",
  "depriv_index_quintile", "period_diag", sep = " + ")

BASE_COVAR_COLS <- c("person_id", "age_cl", "previous_pregnancy", "breast_surgery_2cl",
  "ct_setting_5cl", "rt", "ht_type", "pnuicc_2cl_modified", "her2_status",
  "nb_comor_ql", "Cardiovascular", "Endocrine_and_metabolism",
  "Psychiatric_disorders", "Other", "depriv_index_quintile", "period_diag")

# ── Trial 1 MSM: categorical arm, covariate-adjusted, g-computation ───────────
analyze_survival_outcomes_trial1 <- function(weighted_df, tmax = 72) {

  fit_stratum <- function(dat) {
    dat <- dat %>%
      mutate(time = month, timesq = time^2, outcome = as.numeric(outcome)) %>%
      filter(time <= 72)
    if (!("boot_n" %in% names(dat))) dat$boot_n <- 1

    model_glm <- suppressWarnings(glm(
      as.formula(paste("I(outcome == 0) ~ arm * time + arm * timesq +", BASE_COVARS)),
      data = dat, family = binomial(link = "logit"),
      weights = weight_global_truncated98 * boot_n))

    baseline <- dat %>% filter(time == 0) %>% select(all_of(BASE_COVAR_COLS)) %>% unique()
    regimes  <- unique(dat$arm)

    tidyr::expand_grid(person_id = baseline$person_id, time = 0:tmax, arm = regimes) %>%
      left_join(baseline, by = "person_id") %>%
      mutate(timesq = time^2,
             p_not_event = predict(model_glm, newdata = ., type = "response")) %>%
      arrange(arm, person_id, time) %>%
      group_by(arm, person_id) %>% mutate(p_surv_i = cumprod(p_not_event)) %>% ungroup() %>%
      group_by(arm, time) %>% summarise(p_surv = mean(p_surv_i), .groups = "drop")
  }

  list(predict_df           = fit_stratum(weighted_df),
       predict_df_low_risk  = fit_stratum(weighted_df %>% filter(risk == "Low/Intermediate risk")),
       predict_df_high_risk = fit_stratum(weighted_df %>% filter(risk == "High risk")))
}

# ── Trial 2 MSM: dose-response for the nine windows + control model ───────────
# Interruption arms: interruption timing modelled as a continuous exposure
# (mean ET duration of each window) with linear + quadratic terms interacted
# with follow-up time. Control arms: categorical model. Strategy-specific
# curves are obtained by g-computation and combined.
analyze_survival_outcomes_trial2 <- function(weighted_df, tmax = 72) {

  arm_timing <- data.frame(
    arm = c("x6-12_y3_z1000","x12-18_y3_z1000","x18-24_y3_z1000","x24-30_y3_z1000",
            "x30-36_y3_z1000","x36-42_y3_z1000","x42-48_y3_z1000","x48-54_y3_z1000",
            "x54-60_y3_z1000"),
    interrupt_time = c(9, 15, 21, 27, 33, 39, 45, 51, 57),   # mean ET duration per window
    stringsAsFactors = FALSE)

  fit_stratum <- function(dat) {
    dat <- dat %>%
      mutate(time = month, timesq = time^2, outcome = as.numeric(outcome)) %>%
      filter(time <= 72) %>%
      left_join(arm_timing, by = "arm")
    if (!("boot_n" %in% names(dat))) dat$boot_n <- 1

    # interruption arms: dose-response in interruption timing
    model_int <- suppressWarnings(glm(
      as.formula(paste("I(outcome == 0) ~ (interrupt_time + I(interrupt_time^2)) *",
                       "(time + timesq) +", BASE_COVARS)),
      data = dat %>% filter(is_control == 0), family = binomial(link = "logit"),
      weights = weight_global_truncated98 * boot_n))

    # control arms: categorical
    model_ctrl <- suppressWarnings(glm(
      as.formula(paste("I(outcome == 0) ~ arm * time + arm * timesq +", BASE_COVARS)),
      data = dat %>% filter(is_control == 1), family = binomial(link = "logit"),
      weights = weight_global_truncated98 * boot_n))

    baseline <- dat %>% filter(time == 0) %>% select(all_of(BASE_COVAR_COLS)) %>% unique()

    pred_int <- tidyr::expand_grid(person_id = baseline$person_id, time = 0:tmax,
                                   arm = arm_timing$arm) %>%
      left_join(baseline, by = "person_id") %>%
      left_join(arm_timing, by = "arm") %>%
      mutate(timesq = time^2,
             p_not_event = predict(model_int, newdata = ., type = "response")) %>%
      arrange(arm, person_id, time) %>%
      group_by(arm, person_id) %>% mutate(p_surv_i = cumprod(p_not_event)) %>% ungroup() %>%
      group_by(arm, time) %>% summarise(p_surv = mean(p_surv_i), .groups = "drop")

    pred_ctrl <- tidyr::expand_grid(person_id = baseline$person_id, time = 0:tmax,
                                    arm = c("Control (>= 24mo ET)", "Control (>= 60mo ET)")) %>%
      left_join(baseline, by = "person_id") %>%
      mutate(timesq = time^2,
             p_not_event = predict(model_ctrl, newdata = ., type = "response")) %>%
      arrange(arm, person_id, time) %>%
      group_by(arm, person_id) %>% mutate(p_surv_i = cumprod(p_not_event)) %>% ungroup() %>%
      group_by(arm, time) %>% summarise(p_surv = mean(p_surv_i), .groups = "drop")

    bind_rows(pred_int, pred_ctrl)
  }

  list(predict_df           = fit_stratum(weighted_df),
       predict_df_low_risk  = fit_stratum(weighted_df %>% filter(risk == "Low/Intermediate risk")),
       predict_df_high_risk = fit_stratum(weighted_df %>% filter(risk == "High risk")))
}

# ── bootstrap_survival_analysis ──────────────────────────────────────────────
# Nonparametric bootstrap with stratified cluster resampling: individuals are
# resampled with replacement within each (arm x risk) stratum, preserving the
# original stratum sizes. IP weights are recomputed on each bootstrap sample and
# the trial's MSM (`analyze_fun`) refitted.
bootstrap_survival_analysis <- function(person_month_df, censored_df, analyze_fun,
                                        n_bootstrap = 1000, seed = 123,
                                        id_col = "person_id", arm_col = "arm",
                                        risk_col = "risk", verbose = TRUE) {

  stopifnot(is.function(analyze_fun))
  data.table::setDT(person_month_df)
  data.table::setDT(censored_df)
  set.seed(seed)

  dt <- person_month_df; censor <- censored_df
  person_tbl   <- unique(censor[, c(id_col, arm_col, risk_col), with = FALSE])
  person_tbl   <- person_tbl[!is.na(get(id_col))]
  n_persons    <- nrow(person_tbl)
  strata_sizes <- person_tbl[, .N, by = c(arm_col, risk_col)]

  bootstrap_results <- vector("list", n_bootstrap)
  tic0 <- Sys.time()
  if (verbose) cat("Starting bootstrap with", n_bootstrap, "iterations...\n")

  for (b in seq_len(n_bootstrap)) {
    if (verbose && (b %% 50 == 0 || b == 1))
      cat(sprintf("  bootstrap %d / %d (elapsed %.1f min)\n", b, n_bootstrap,
                  as.numeric(difftime(Sys.time(), tic0, units = "mins"))))

    bootstrap_results[[b]] <- suppressMessages(tryCatch({

      sampled_list <- vector("list", nrow(strata_sizes))
      for (k in seq_len(nrow(strata_sizes))) {
        a_val <- strata_sizes[[arm_col]][k]; r_val <- strata_sizes[[risk_col]][k]
        ids_k <- person_tbl[get(arm_col) == a_val & get(risk_col) == r_val, get(id_col)]
        sampled_list[[k]] <- sample(ids_k, size = strata_sizes[["N"]][k], replace = TRUE)
      }
      sampled_ids <- unlist(sampled_list, use.names = FALSE)

      boot_id_dt <- data.table::as.data.table(table(sampled_ids))
      data.table::setnames(boot_id_dt, c(id_col, "boot_n"))
      boot_id_dt[, boot_n := as.integer(boot_n)]

      boot_dt   <- boot_id_dt[dt,     on = id_col, nomatch = 0L]
      boot_cens <- boot_id_dt[censor, on = id_col, nomatch = 0L]

      weighted <- calculate_ip_weights(person_month_df = boot_dt, censored_df = boot_cens)
      if (!("boot_n" %in% names(weighted))) weighted[, boot_n := boot_dt$boot_n]

      surv_out <- analyze_fun(weighted)

      list(iteration = b, convergence = TRUE,
           predict_df = surv_out$predict_df,
           predict_df_low_risk = surv_out$predict_df_low_risk,
           predict_df_high_risk = surv_out$predict_df_high_risk, error = NULL)

    }, error = function(e) {
      list(iteration = b, convergence = FALSE, predict_df = NULL,
           predict_df_low_risk = NULL, predict_df_high_risk = NULL,
           error = conditionMessage(e))
    }))

    rm(list = intersect(ls(), c("boot_dt","boot_cens","weighted","surv_out",
                                "boot_id_dt","sampled_list","sampled_ids")))
    gc(FALSE)
  }

  n_success <- sum(vapply(bootstrap_results, `[[`, logical(1), "convergence"))
  if (verbose) cat("Completed:", n_success, "/", n_bootstrap, "successful.\n")

  structure(list(results = bootstrap_results, n_bootstrap = n_bootstrap,
                 n_successful = n_success, original_n_persons = n_persons,
                 settings = list(seed = seed)),
            class = "bootstrap_survival")
}

# ── extract_ci_for_plot ──────────────────────────────────────────────────────
# Bootstrap percentile CIs for the cumulative-incidence curves (for figures).
extract_ci_for_plot <- function(bootstrap_results, confidence_level = 0.95,
                                prediction_col = "p_surv",
                                which_dfs = c("predict_df", "predict_df_low_risk",
                                              "predict_df_high_risk"),
                                scale_to_pct = FALSE) {

  ok <- purrr::map_lgl(bootstrap_results$results, ~ isTRUE(.x$convergence))
  successful <- bootstrap_results$results[ok]

  pull_pred <- function(slot) {
    preds <- purrr::map(successful, ~ .x[[slot]])
    preds <- preds[!purrr::map_lgl(preds, is.null)]
    if (!length(preds)) return(NULL)
    bind_rows(preds, .id = "bootstrap_id") %>% mutate(pred_set = slot)
  }

  all_pred <- purrr::map(which_dfs, pull_pred) %>% purrr::compact() %>% bind_rows()
  if (!nrow(all_pred)) stop("No prediction data frames found.")
  if (scale_to_pct) all_pred <- all_pred %>% mutate(!!prediction_col := .data[[prediction_col]] * 100)

  a <- 1 - confidence_level
  all_pred %>%
    group_by(time, arm, pred_set) %>%
    summarise(n_bootstrap = n(),
              mean_pred = mean(.data[[prediction_col]], na.rm = TRUE),
              median_pred = median(.data[[prediction_col]], na.rm = TRUE),
              ci_lower = quantile(.data[[prediction_col]], a/2, na.rm = TRUE),
              ci_upper = quantile(.data[[prediction_col]], 1 - a/2, na.rm = TRUE),
              .groups = "drop")
}

# ── extract_5year_difference ─────────────────────────────────────────────────
# Per-arm 5-year cumulative incidence and 5-year risk difference vs each control,
# with bootstrap percentile CIs. diff_pct = control - interruption (positive =
# higher risk with interruption). CI text uses the "xx to xx" format.
extract_5year_difference <- function(bootstrap_results, confidence_level = 0.95,
                                     target_time = 60, control_pattern = "^Control") {

  successful <- bootstrap_results$results[
    purrr::map_lgl(bootstrap_results$results, ~ isTRUE(.x$convergence))]
  if (!length(successful)) stop("No successful bootstrap iterations.")

  extract_5y <- function(df_pred) {
    df_pred %>% group_by(arm) %>%
      mutate(time_diff = abs(time - target_time)) %>%
      filter(time_diff == min(time_diff)) %>% select(arm, time, p_surv) %>% ungroup()
  }
  a <- 1 - confidence_level; lower_q <- a/2; upper_q <- 1 - a/2

  compute_for_dataset <- function(predictions, dataset_label) {
    predictions <- predictions[!purrr::map_lgl(predictions, is.null)]
    if (!length(predictions)) return(NULL)
    all_5y <- purrr::map_dfr(predictions, extract_5y, .id = "bootstrap_id")

    per_arm <- all_5y %>% group_by(arm) %>%
      summarise(n_bootstrap = n(), actual_time = round(mean(time), 2),
                median_risk_pct = round(median((1 - p_surv) * 100), 1),
                risk_lower_pct  = round(quantile((1 - p_surv) * 100, lower_q), 1),
                risk_upper_pct  = round(quantile((1 - p_surv) * 100, upper_q), 1),
                .groups = "drop") %>%
      mutate(risk_formatted = sprintf("%.1f (%.1f to %.1f)",
                                      median_risk_pct, risk_lower_pct, risk_upper_pct),
             dataset = dataset_label)

    control_arms <- unique(all_5y$arm[grepl(control_pattern, all_5y$arm)])
    non_control  <- setdiff(unique(all_5y$arm), control_arms)
    base_dat <- all_5y %>% select(bootstrap_id, arm, p_surv)

    diff_results <- base_dat %>% filter(arm %in% non_control) %>%
      rename(test_arm = arm, p_surv_test = p_surv) %>%
      inner_join(base_dat %>% filter(arm %in% control_arms) %>%
                   rename(control_arm = arm, p_surv_control = p_surv), by = "bootstrap_id") %>%
      mutate(diff_pct = (p_surv_control - p_surv_test) * 100) %>%
      group_by(test_arm, control_arm) %>%
      summarise(n_bootstrap = n(),
                median_diff_pct = median(diff_pct, na.rm = TRUE),
                diff_lower_pct  = quantile(diff_pct, lower_q, na.rm = TRUE),
                diff_upper_pct  = quantile(diff_pct, upper_q, na.rm = TRUE), .groups = "drop") %>%
      mutate(diff_formatted = sprintf("%.1f (%.1f to %.1f)",
                                      median_diff_pct, diff_lower_pct, diff_upper_pct),
             dataset = dataset_label)

    list(per_arm = per_arm, diff_vs_control = diff_results)
  }

  dataset_map <- c(overall = "predict_df", low_risk = "predict_df_low_risk",
                   high_risk = "predict_df_high_risk")
  res_list <- purrr::imap(dataset_map, function(df_name, label)
    compute_for_dataset(purrr::map(successful, ~ .x[[df_name]]), label)) %>% purrr::compact()

  list(per_arm = purrr::map_dfr(res_list, "per_arm"),
       diff_vs_control = purrr::map_dfr(res_list, "diff_vs_control"))
}
