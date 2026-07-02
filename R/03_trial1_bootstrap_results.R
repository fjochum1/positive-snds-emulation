################################################################################
# 03_trial1_bootstrap_results.R
#
# TRIAL 1 — bootstrap confidence intervals and 5-year risk differences.
# Nonparametric bootstrap (1000 replicates): IP weights recomputed and MSM 
# refitted on each replicate.
# Outputs: per-arm 5-year risks and interruption-vs-control 5-year risk
# differences, with 95% percentile CIs
#
# Run 01_trial1_clone_censor_weight.R first (defines calculate_ip_weights and
# builds censored_df).
################################################################################

source("R/00_setup.R")
source("R/functions_msm_bootstrap.R")

################################################################################
# Bootstrap
################################################################################

set.seed(123)  # reproducibility
bootstrap_first <- bootstrap_survival_analysis(
  person_month_df = pop_positive_first,
  censored_df     = censored_df,
  analyze_fun     = analyze_survival_outcomes_trial1,   # adjusted MSM, g-computation
  n_bootstrap     = 1000,
  seed            = 123)

save(bootstrap_first, file = file.path(results_dir, "bootstrap_first.RData"))

################################################################################
# 5-year risks and risk differences (5-year horizon = month 60 from time zero)
################################################################################

res     <- extract_5year_difference(bootstrap_first, target_time = 60)
df_arm  <- res$per_arm
df_diff <- res$diff_vs_control

# Readable labels
recode_dataset <- function(x) factor(x, levels = c("overall", "low_risk", "high_risk"),
                                      labels = c("Overall", "Low/Intermediate risk", "High risk"))
recode_control <- function(x) factor(x, levels = c("Control (>= 24mo ET)",),
                                      labels = c("Control >= 60mo"))

df_diff <- df_diff %>%
  mutate(interruption_window = "18-30mo",             # single window in Trial 1
         dataset = recode_dataset(dataset),
         control_arm = recode_control(control_arm)) %>%
  arrange(control_arm, dataset) %>%
  select(dataset, control_arm, interruption_window,
         median_diff_pct, diff_lower_pct, diff_upper_pct, diff_formatted, n_bootstrap)

df_arm <- df_arm %>% mutate(dataset = recode_dataset(dataset)) %>% arrange(dataset, arm)

# Save results 
write.csv(df_diff, file.path(results_dir, "trial1_5y_risk_differences.csv"), row.names = FALSE)
write.csv(df_arm,  file.path(results_dir, "trial1_5y_risks_by_arm.csv"),      row.names = FALSE)

print(df_diff)
cat("\nTrial 1 bootstrap results saved to ", results_dir, "\n")


