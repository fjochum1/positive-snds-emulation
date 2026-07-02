################################################################################
# 06_sensitivity_analyses.R
#
# Sensitivity analyses (Trial 2 unless stated). Each is a SINGLE modification to
# the main pipeline; everything else (clone -> censor -> IP weights -> MSM ->
# bootstrap -> 5-year risk differences) is unchanged. For each analysis only the
# change is shown, then the standard pipeline is re-run:
#
#   cloned_df   <- clone_for_interventions(pop, interventions)
#   censored_df <- apply_censoring_rules(cloned_df)
#   weighted_df <- calculate_ip_weights(pop, censored_df)
#   boot        <- bootstrap_survival_analysis(pop, censored_df,
#                                              analyze_fun = analyze_survival_outcomes_trial2)
#   res         <- extract_5year_difference(boot, target_time = 60)
#
################################################################################

source("R/00_setup.R")
source("R/functions_msm_bootstrap.R")

# ==============================================================================
# A. Pregnancy-seeking population
#    CHANGE: additionally censor women with no pregnancy-seeking indicator in the
#    window [interruption - 6 mo, interruption + 24 mo], and add a pregnancy-
#    seeking censoring weight (extra IPCW factor). Re-run the pipeline.
# ------------------------------------------------------------------------------
#   In apply_censoring_rules(), add:
#     censor_no_pregnancy_seeking <- interventions_df %>%
#       filter(month == et_stop_month + 24 & planning_pregnancy_window == 0) %>%
#       transmute(person_id, arm, censor_month = month,
#                 censor_reason = "no_pregnancy_seeking")
#   and include it in bind_rows(...). In calculate_ip_weights(), multiply
#   ipcw_weight by an extra factor from a pooled-logistic model of pregnancy-
#   seeking (same covariates), analogous to weight1/weight4.

# ==============================================================================
# B. Enforced ET resumption (within the pregnancy-seeking population)
#    CHANGE: on top of (A), censor women who did not resume ET within 24 months
#    after interruption (or 6 months postpartum); reuse the ET-utilisation model
#    (Pr2) for the extra censoring weight.
# ------------------------------------------------------------------------------
#   In apply_censoring_rules(), add for the interruption arms:
#     censor_no_resumption <- interventions_df %>%
#       filter(month == pmin(et_stop_month + 24, pregnancy_end_month + 6) & et == 0) %>%
#       transmute(person_id, arm, censor_month = month,
#                 censor_reason = "no_resumption")
#   and a corresponding 1/Pr2 weight factor for the resumption window.

# ==============================================================================
# C. Negative-control outcome (otolaryngologist visit)
#    CHANGE: swap the outcome (breast-cancer event) for the negative-control
#    event; keep the same clones, censoring and weights.
# ------------------------------------------------------------------------------
#   In apply_censoring_rules(), replace:
#     mutate(outcome = relapse, outcome_month = relapse_month)
#   by:
#     mutate(outcome = orl_visit, outcome_month = orl_visit_month)
#   Then re-run calculate_ip_weights() -> bootstrap -> extract_5year_difference().

# ==============================================================================
# D. Alternative weight truncation (95th instead of 98th percentile)
#    CHANGE: use the 95th-percentile-truncated weight in the MSM.
# ------------------------------------------------------------------------------
#   In analyze_survival_outcomes_trial2() (and _trial1), replace every
#     weights = weight_global_truncated98 * boot_n
#   by
#     weights = weight_global_truncated95 * boot_n

# ==============================================================================
# E. Alternative dose-response specification (restricted cubic splines)
#    CHANGE: model interruption timing with natural splines instead of a
#    linear + quadratic polynomial, in the Trial 2 intervention MSM.
# ------------------------------------------------------------------------------
#   In analyze_survival_outcomes_trial2(), replace the intervention formula
#     I(outcome == 0) ~ (interrupt_time + I(interrupt_time^2)) * (time + timesq) + <covariates>
#   by
#     I(outcome == 0) ~ ns(interrupt_time, df = 3) * (time + timesq) + <covariates>
#   (knots at quantiles; all covariates, weights, strata and controls unchanged).

# ==============================================================================
# F. Alternative risk classification: nodal status alone
#    CHANGE: redefine the stratifying variable `risk` as node-positive vs
#    node-negative, then re-run (strata are read from `risk`).
# ------------------------------------------------------------------------------
#   pop_positive_second <- pop_positive_second %>%
#     mutate(risk = if_else(pnuicc_2cl_modified == "Node-positive",
#                           "High risk", "Low/Intermediate risk"))
#   # then clone -> censor -> weight -> bootstrap -> extract, as usual.

# ==============================================================================
# G. Age-stratified analysis (<35 vs >=35 years)
#    CHANGE: run the whole pipeline within each age subgroup (the "Overall"
#    output of each run gives that age group's estimate).
# ------------------------------------------------------------------------------
#   for (grp in list(`<35` = quote(age < 35), `>=35` = quote(age >= 35))) {
#     pop_g  <- pop_positive_second %>% filter(!!grp)
#     cens_g <- apply_censoring_rules(clone_for_interventions(pop_g, interventions))
#     boot_g <- bootstrap_survival_analysis(pop_g, cens_g,
#                                           analyze_fun = analyze_survival_outcomes_trial2)
#     # extract_5year_difference(boot_g)$diff_vs_control  -> that age group
#   }

# ==============================================================================
# H. Dynamic continuation strategy (adverse-event-related interruption allowed):
#    in the CONTROL arms only, once a serious ET-related adverse event
#    has occurred the woman is released from the continuation constraint (she may
#    stop or continue ET thereafter) -- her discontinuation is no longer a
#    deviation. The interruption arms are unaffected (an early stop there defines
#    an earlier window, not a deviation).
# ------------------------------------------------------------------------------
#
#   (i)  Censoring -- in apply_censoring_rules(), for the control ET-stop rules
#        (censor_et_stop_before_m24 / _m60) add:  filter(ae_ever == 0)
#        so an AE-driven discontinuation is not treated as a deviation.
#
#   (ii) Weights -- in calculate_ip_weights(), add a case to weight5 (before the
#        et cases) so a released control clone is no longer at risk (weight 1):
#          weight5 = case_when(
#            month >= x_min ~ 1,
#            is_control == 1 & ae_ever == 1 ~ 1,   # released after AE -> weight 1
#            month < x_min & et == 0 ~ 0,
#            month < x_min & et == 1 ~ 1 / Pr2,
#            TRUE ~ NA_real_)
#


