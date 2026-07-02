################################################################################
# 02_trial2_clone_censor_weight.R
#
# TRIAL 2 — extension to the timing of ET interruption.
# Time zero = ET initiation.
# Nine interruption windows (6-12, 12-18, ..., 54-60 months of ET) vs two
# continuation controls (uninterrupted ET for >= 24 or >= 60 months).
# Pipeline: load person-month data -> clone -> censor at first deviation ->
#           inverse-probability-of-censoring weights -> save weighted_df.
#
################################################################################

source("R/00_setup.R")

################################################################################
# Load data
################################################################################

load(file.path(data_dir, "pop_positive_second.RData"))  # person-month data, time zero = ET initiation

cat("Unique individuals:", length(unique(pop_positive_second$person_id)), "\n") # 22 916 women

################################################################################
# Definition of the interventions (nine 6-month interruption windows)
################################################################################

create_intervention_grid <- function() {
  x_pairs <- data.frame(
    x_min = c(6, 12, 18, 24, 30, 36, 42, 48, 54),
    x_max = c(12, 18, 24, 30, 36, 42, 48, 54, 60))
  x_pairs %>% crossing(y = c(3), z = c(1000))
}

interventions <- create_intervention_grid()
cat("Total number of intervention arms:", nrow(interventions), "\n")

################################################################################
# CLONE
################################################################################
# Each individual is cloned to the nine interruption windows plus two
# continuation controls (uninterrupted ET for >= 24 and >= 60 months).

clone_for_interventions <- function(data, interventions) {

  control_data1 <- data %>%
    mutate(arm_id = "Control (>= 24mo ET)", x_min = 24, x_max = NA, y = NA, z = NA,
           is_control = 1, arm = "Control (>= 24mo ET)")

  control_data2 <- data %>%
    mutate(arm_id = "Control (>= 60mo ET)", x_min = 60, x_max = NA, y = NA, z = NA,
           is_control = 1, arm = "Control (>= 60mo ET)")

  intervention_data <- NULL
  for (i in seq_len(nrow(interventions))) {
    int <- interventions[i, ]
    clone <- data %>%
      mutate(arm_id = paste0("arm_", i),
             x_min = int$x_min, x_max = int$x_max, y = int$y, z = int$z,
             is_control = 0,
             arm = paste0("x", x_min, "-", x_max, "_y", y, "_z", z))
    intervention_data <- bind_rows(intervention_data, clone)
  }

  bind_rows(control_data1, control_data2, intervention_data)
}

cloned_df <- clone_for_interventions(pop_positive_second, interventions)

################################################################################
# CENSOR
################################################################################
# Interruption arms: ET stopped before x_min; not stopped by x_max; ET resumed
# during the 3-month washout; pregnancy before the end of washout.
# Control arms: ET stopped before the required horizon (24 or 60 months); or
# pregnancy before 60 months of ET.

apply_censoring_rules <- function(cloned_df) {

  df <- cloned_df
  interventions_df <- df %>% filter(is_control == 0)

  et_stop <- interventions_df %>%
    filter(et == 0) %>%
    group_by(person_id) %>% filter(month == min(month)) %>% ungroup() %>%
    rename(et_stop_month = month) %>% select(person_id, et_stop_month)

  interventions_df <- interventions_df %>% left_join(et_stop, by = "person_id")

  censor_ET_stop_before_xmin <- interventions_df %>%
    filter(et == 0 & month < x_min) %>%
    group_by(person_id, arm) %>% filter(month == min(month)) %>% ungroup() %>%
    mutate(censor_month = month, censor_reason = "ET_stop_before_xmin") %>%
    select(person_id, arm, censor_month, censor_reason) %>% unique()

  censor_at_xmax <- interventions_df %>%
    filter(is.na(et_stop_month) | et_stop_month >= x_max) %>%
    mutate(censor_month = x_max, censor_reason = "no_ET_stop_by_xmax") %>%
    select(person_id, arm, censor_month, censor_reason) %>% unique()

  censor_restart_early <- interventions_df %>%
    filter(!is.na(et_stop_month)) %>%
    filter((et == 1 & month == et_stop_month + 1) |
           (et == 1 & et_lag1 == 0 & month == et_stop_month + 2)) %>%
    mutate(censor_month = month, censor_reason = "et_restart_early") %>%
    select(person_id, arm, censor_month, censor_reason) %>% unique()

  censor_pregnancy_without_washout <- interventions_df %>%
    filter(pregnancy == 1 & month < (et_stop_month + 3)) %>%
    select(person_id, month, arm) %>% unique() %>%
    group_by(person_id, arm) %>% filter(month == min(month)) %>% ungroup() %>%
    rename(censor_month = month) %>% mutate(censor_reason = "pregnancy_without_washout")

  # (control) ET stopped before required horizon
  censor_et_stop_before_m24 <- df %>%
    filter(arm == "Control (>= 24mo ET)") %>% filter(et == 0 & month < 24) %>%
    select(person_id, month, arm) %>% unique() %>%
    group_by(person_id) %>% filter(month == min(month)) %>% ungroup() %>%
    rename(censor_month = month) %>% mutate(censor_reason = "et_stop_before_m24")

  censor_et_stop_before_m60 <- df %>%
    filter(arm == "Control (>= 60mo ET)") %>% filter(et == 0 & month < 60) %>%
    select(person_id, month, arm) %>% unique() %>%
    group_by(person_id) %>% filter(month == min(month)) %>% ungroup() %>%
    rename(censor_month = month) %>% mutate(censor_reason = "et_stop_before_m60")

  # (control) pregnancy before 60 months of ET
  censor_pregnancy_before_m60 <- df %>%
    filter(is_control == 1) %>% filter(pregnancy == 1 & month < 60) %>%
    select(person_id, month, arm) %>% unique() %>%
    group_by(person_id) %>% filter(month == min(month)) %>% ungroup() %>%
    rename(censor_month = month) %>% mutate(censor_reason = "pregnancy_before_m60")

  all_censoring <- bind_rows(
      censor_ET_stop_before_xmin, censor_at_xmax, censor_restart_early,
      censor_pregnancy_without_washout, censor_et_stop_before_m24,
      censor_et_stop_before_m60, censor_pregnancy_before_m60) %>%
    group_by(person_id, arm) %>%
    filter(censor_month == min(censor_month, na.rm = TRUE)) %>% unique() %>%
    mutate(censor_priority = case_when(
      censor_reason == "pregnancy_without_washout" ~ 1,
      censor_reason == "no_ET_stop_by_xmax"        ~ 2,
      TRUE                                         ~ 3)) %>%
    group_by(person_id, arm, censor_month) %>% arrange(censor_priority) %>%
    slice(1) %>% ungroup() %>% select(-censor_priority)

  censored_df <- df %>%
    left_join(all_censoring, by = c("person_id", "arm")) %>%
    mutate(outcome = relapse, outcome_month = relapse_month) %>%
    mutate(outcome = ifelse(
      (!is.na(outcome_month) & outcome_month <= censor_month & outcome == 1) |
      (is.na(censor_month) & !is.na(outcome_month) & outcome == 1), 1, 0)) %>%
    mutate(censored = !is.na(censor_month) & month >= censor_month,
           censored = ifelse(!is.na(outcome_month) & outcome_month <= censor_month & outcome == 1,
                             FALSE, censored),
           censor_reason = if_else(!censored, NA, censor_reason),
           censor_month  = if_else(!censored, NA, censor_month)) %>%
    group_by(person_id, arm) %>%
    filter((!is.na(censor_month) & month <= censor_month) | is.na(censor_month)) %>%
    filter((!is.na(outcome_month) & month <= outcome_month) | is.na(outcome_month)) %>%
    ungroup()

  variable_to_keep <- interventions_df %>%
    select(person_id, arm, et_stop_month) %>% unique()

  left_join(censored_df, variable_to_keep)
}

create_censoring_summary <- function(censored_df) {
  censored_types <- as.data.frame(censored_df) %>%
    group_by(person_id, arm) %>%
    summarise(
      arm = first(arm),
      censor_month = if (all(is.na(censor_month))) NA else first(censor_month[!is.na(censor_month)]),
      last_month = max(month),
      censored = any(censored),
      censor_ET_stop_before_xmin       = any(censor_reason == "ET_stop_before_xmin" & censored),
      censor_pregnancy_without_washout = any(censor_reason == "pregnancy_without_washout" & censored),
      censor_no_ET_stop_by_xmax        = any(censor_reason == "no_ET_stop_by_xmax" & censored),
      censor_et_restart_early          = any(censor_reason == "et_restart_early" & censored),
      censor_et_stop_before_m24        = any(censor_reason == "et_stop_before_m24" & censored),
      censor_et_stop_before_m60        = any(censor_reason == "et_stop_before_m60" & censored),
      censor_pregnancy_before_m60      = any(censor_reason == "pregnancy_before_m60" & censored),
      outcome = any(outcome == 1), .groups = "drop") %>%
    mutate(admin_censored = if_else(!censored & outcome == 0, TRUE, FALSE))

  censored_types %>%
    group_by(arm) %>%
    summarise(
      total_patients = n(),
      n_censored_ET_stop_before_xmin       = sum(censor_ET_stop_before_xmin, na.rm = TRUE),
      n_censored_pregnancy_without_washout = sum(censor_pregnancy_without_washout, na.rm = TRUE),
      n_censored_no_ET_stop_by_xmax        = sum(censor_no_ET_stop_by_xmax, na.rm = TRUE),
      n_censored_et_restart_early          = sum(censor_et_restart_early, na.rm = TRUE),
      n_censored_et_stop_before_m24        = sum(censor_et_stop_before_m24, na.rm = TRUE),
      n_censored_et_stop_before_m60        = sum(censor_et_stop_before_m60, na.rm = TRUE),
      n_censored_pregnancy_before_m60      = sum(censor_pregnancy_before_m60, na.rm = TRUE),
      n_admin_censored = sum(admin_censored, na.rm = TRUE),
      n_outcome = sum(outcome), .groups = "drop") %>%
    arrange(arm)
}

censored_df   <- apply_censoring_rules(cloned_df)
summary_table <- create_censoring_summary(censored_df)
print(summary_table)
save(censored_df, file = file.path(results_dir, "censored_df_second.RData"))

################################################################################
# WEIGHT (inverse-probability-of-censoring weights)
################################################################################
#   Model 1 (Pr1): monthly probability of pregnancy
#   Model 2 (Pr2): monthly probability of being on ET

calculate_ip_weights <- function(person_month_df, censored_df) {

  df <- person_month_df %>%
    mutate(time = month, timesq = time^2,
           time_cat = cut(time, breaks = c(0,6,18,30,42,60,120),
                          include.lowest = TRUE, right = TRUE))

  categorical_vars <- c(
    "age_cl", "previous_pregnancy", "breast_surgery_2cl", "pnuicc_2cl_modified",
    "ct_setting_5cl", "rt", "ht_type", "her2_status", "nb_comor_ql",
    "Cardiovascular", "Endocrine_and_metabolism", "Psychiatric_disorders", "Other",
    "depriv_index_quintile", "planning_pregnancy_past_6mo",
    "et", "et_lag1", "et_lag2", "et_lag3", "followup_past6m", "time_cat",
    "pregnancy", "had_pregnancy", "period_diag")

  df <- df %>% mutate(across(all_of(categorical_vars), as.character))

  first_pregnancy_df <- df %>%
    filter(pregnancy == 1) %>%
    group_by(person_id) %>% summarise(first_pregnancy_month = min(month), .groups = "drop")

  df <- df %>% left_join(first_pregnancy_df, by = "person_id")

  model1_data <- df %>%
    filter(is.na(first_pregnancy_month) | month <= first_pregnancy_month) %>%
    mutate(outcome_pregnancy = as.integer(month == first_pregnancy_month),
           outcome_pregnancy = if_else(is.na(outcome_pregnancy), 0L, outcome_pregnancy)) %>%
    select(person_id, month, outcome_pregnancy, all_of(categorical_vars),
           cumulative_side_effects) %>% unique()

  base_vars <- c("age_cl", "previous_pregnancy", "breast_surgery_2cl",
                 "pnuicc_2cl_modified", "ct_setting_5cl", "rt", "ht_type",
                 "her2_status", "nb_comor_ql", "Cardiovascular",
                 "Endocrine_and_metabolism", "Psychiatric_disorders", "Other",
                 "depriv_index_quintile", "period_diag")
  base_pred <- paste(base_vars, collapse = " + ")

  # Model 1: probability of pregnancy (Pr1)
  var_model1 <- paste0(base_pred, " + ns(month, df=2) + cumulative_side_effects",
                       " + followup_past6m + et + et_lag1 + et_lag2 + planning_pregnancy_past_6mo")
  # Model 2: probability of being on ET (Pr2)
  var_model2 <- paste0(base_pred, " + ns(month, df=2) + et_lag1 + et_lag2",
                       " + cumulative_side_effects + followup_past6m",
                       " + planning_pregnancy_past_6mo + pregnancy + had_pregnancy")

  fit <- list(
    model1 = speedglm(as.formula(paste("I(outcome_pregnancy == 1) ~", var_model1)),
                      data = model1_data, family = binomial()),
    model2 = speedglm(as.formula(paste("I(et == 1) ~", var_model2)),
                      data = df, family = binomial()))

  tidy(fit$model1, exponentiate = TRUE, conf.int = TRUE) %>%
    select(-std.error, -p.value, -statistic) %>% print(n = Inf)
  tidy(fit$model2, exponentiate = TRUE, conf.int = TRUE) %>%
    select(-std.error, -p.value, -statistic) %>% print(n = Inf)

  Prob <- df %>%
    mutate(Pr1 = predict(fit$model1, df, type = "response"),
           Pr2 = predict(fit$model2, df, type = "response")) %>%
    select(person_id, month, Pr1, Pr2)

  df <- censored_df %>%
    left_join(first_pregnancy_df, by = "person_id") %>%
    left_join(Prob)

  df <- df %>%
    mutate(
      # pregnancy before end of washout (interruption arms)
      weight1 = case_when(
        is_control == 1 ~ NA_real_,
        month >= et_stop_month + 3 ~ 1,
        pregnancy == 1 & month < et_stop_month + 3 ~ 0,
        pregnancy == 0 & month < et_stop_month + 3 ~ 1 / (1 - Pr1),
        TRUE ~ NA_real_),
      # not stopped at x_max (interruption arms)
      weight2 = case_when(
        is_control == 1 ~ NA_real_,
        month != x_max ~ 1,
        month == x_max & et_stop_month < x_max ~ 1 / (1 - Pr2),
        month == x_max & (et_stop_month >= x_max | is.na(et_stop_month)) ~ 0,
        TRUE ~ NA_real_),
      # ET resumed during washout (interruption arms)
      weight3 = case_when(
        is_control == 1 ~ NA_real_,
        month == et_stop_month + 1 & et == 0 ~ 1 / (1 - Pr2),
        month == et_stop_month + 2 & et_lag1 == 0 & et == 0 ~ 1 / (1 - Pr2),
        month == et_stop_month + 1 & et == 1 ~ 0,
        month == et_stop_month + 2 & et_lag1 == 0 & et == 1 ~ 0,
        TRUE ~ 1),
      # pregnancy before 60 months of ET (control arms)
      weight4 = case_when(
        is_control == 0 ~ NA_real_,
        month >= 60 ~ 1,
        month < 60 & month == first_pregnancy_month ~ 0,
        month < 60 & is.na(first_pregnancy_month) ~ 1 / (1 - Pr1),
        month < 60 & month < first_pregnancy_month ~ 1 / (1 - Pr1),
        TRUE ~ NA_real_),
      # ET stopped before x_min
      weight5 = case_when(
        month >= x_min ~ 1,
        month < x_min & et == 0 ~ 0,
        month < x_min & et == 1 ~ 1 / Pr2,
        TRUE ~ NA_real_),
      ipcw_weight = case_when(
        is_control == 0 ~ weight1 * weight2 * weight3 * weight5,
        is_control == 1 ~ weight4 * weight5)) %>%
    arrange(person_id, arm, month) %>%
    group_by(person_id, arm) %>% mutate(cum_ipcw = cumprod(ipcw_weight)) %>% ungroup()

  out <- df %>%
    group_by(arm) %>%
    mutate(weight_global_truncated98 = pmin(cum_ipcw, quantile(cum_ipcw, 0.98, na.rm = TRUE))) %>%
    ungroup()

  cat("\nCumulative weight (untruncated), max per arm:\n")
  print(aggregate(cum_ipcw ~ arm, data = out, max))
  cat("\nTruncated (98th pct) weights, summary:\n")
  print(summary(out$weight_global_truncated98))

  out
}

weighted_df <- calculate_ip_weights(pop_positive_second, censored_df)
save(weighted_df, file = file.path(results_dir, "weighted_df_second.RData"))

cat("\nTrial 2 clone-censor-weight complete. weighted_df saved to ", results_dir, "\n")
