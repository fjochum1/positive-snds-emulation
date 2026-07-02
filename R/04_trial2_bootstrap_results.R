################################################################################
# 04_trial2_bootstrap_results.R
#
# TRIAL 2 — bootstrap confidence intervals and 5-year risk differences across
# the nine interruption-timing windows, vs each continuation control (>= 24 and
# >= 60 months). Nonparametric bootstrap (1000 replicates): 
# IP weights recomputed and MSM refitted on each replicate.
# Outputs: per-arm 5-year risks and per-window 5-year risk differences with 95%
# percentile CIs 
#
# Run 02_trial2_clone_censor_weight.R first (defines calculate_ip_weights and
# builds censored_df).
################################################################################

source("R/00_setup.R")
source("R/functions_msm_bootstrap.R")

################################################################################
# Bootstrap
################################################################################

set.seed(123)
bootstrap_second <- bootstrap_survival_analysis(
  person_month_df = pop_positive_second,
  censored_df     = censored_df,
  analyze_fun     = analyze_survival_outcomes_trial2,   # dose-response + control MSM, g-computation
  n_bootstrap     = 1000,
  seed            = 123)

save(bootstrap_second, file = file.path(results_dir, "bootstrap_second.RData"))

################################################################################
# 5-year risks and risk differences (5-year horizon = month 60)
################################################################################

res     <- extract_5year_difference(bootstrap_second, target_time = 60)
df_arm  <- res$per_arm
df_diff <- res$diff_vs_control

recode_dataset <- function(x) factor(x, levels = c("overall", "low_risk", "high_risk"),
                                      labels = c("Overall", "Low/Intermediate risk", "High risk"))
recode_control <- function(x) factor(x, levels = c("Control (>= 24mo ET)", "Control (>= 60mo ET)"),
                                      labels = c("Control >= 24mo", "Control >= 60mo"))

# Interruption window label, e.g. "x6-12_y3_z1000" -> "6-12"
window_order <- c("6-12","12-18","18-24","24-30","30-36","36-42","42-48","48-54","54-60")

df_diff <- df_diff %>%
  mutate(interruption_window = factor(gsub("^x([0-9]+-[0-9]+)_.*$", "\\1", test_arm),
                                      levels = window_order),
         dataset = recode_dataset(dataset),
         control_arm = recode_control(control_arm)) %>%
  arrange(control_arm, dataset, interruption_window) %>%
  select(dataset, control_arm, interruption_window,
         median_diff_pct, diff_lower_pct, diff_upper_pct, diff_formatted, n_bootstrap)

df_arm <- df_arm %>% mutate(dataset = recode_dataset(dataset)) %>% arrange(dataset, arm)

write.csv(df_diff, file.path(results_dir, "trial2_5y_risk_differences.csv"), row.names = FALSE)
write.csv(df_arm,  file.path(results_dir, "trial2_5y_risks_by_arm.csv"),      row.names = FALSE)

print(df_diff, n = Inf)
cat("\nTrial 2 bootstrap results saved to ", results_dir, "\n")

