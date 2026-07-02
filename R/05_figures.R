################################################################################
# 05_figures.R
#
# Manuscript figures, regenerated from the bootstrap outputs:
#   Figure 2 — Trial 1 weighted cumulative-incidence curves (interruption vs
#              uninterrupted ET >= 24 months), overall and by baseline risk.
#   Figure 3 — Trial 2 weighted cumulative-incidence curves across the nine
#              interruption windows, overall and by baseline risk.
#   Figure 4 — Trial 2 five-year absolute risk differences (forest plots) vs
#              uninterrupted ET >= 24 and >= 60 months.
#
# Inputs: bootstrap_first / bootstrap_second (from 03_/04_).
# The pregnancy-seeking panels (Figures 2B, 3B, 4C-D) are produced by re-running
# the same functions on the pregnancy-seeking bootstrap objects.
#
# These scripts reproduce the content of the published figures from the tidy
# bootstrap results; exact fonts/colours can be adjusted to taste.
################################################################################

source("R/00_setup.R")
source("R/functions_msm_bootstrap.R")

load(file.path(results_dir, "bootstrap_first.RData"))   # bootstrap_first
load(file.path(results_dir, "bootstrap_second.RData"))  # bootstrap_second

stratum_labels <- c(predict_df = "Overall",
                    predict_df_low_risk = "Low/Intermediate risk",
                    predict_df_high_risk = "High risk")

# Red gradient for interruption windows (earliest = darkest), grey for controls
red_gradient <- function(n) grDevices::colorRampPalette(c("#67000D", "#EF3B2C", "#FCAE91"))(n)

# Turn an arm id into a short label ("x6-12_y3_z1000" -> "6-12"; controls kept)
relabel_arm <- function(arm, trial1 = FALSE) {
  out <- ifelse(grepl("^Control .*24", arm), "Control ≥ 24mo",
         ifelse(grepl("^Control .*60", arm), "Control ≥ 60mo",
                sub("^x([0-9]+-[0-9]+)_.*$", "\\1", arm)))
  if (trial1) out[out == "0-12"] <- "18-30mo"
  out
}

################################################################################
# Cumulative-incidence curves (Figures 2 and 3)
################################################################################

plot_cuminc <- function(boot, trial1 = FALSE, tmax = 72,
                        keep_controls = c("Control ≥ 24mo", "Control ≥ 60mo")) {

  ci <- extract_ci_for_plot(boot, which_dfs = names(stratum_labels),
                            scale_to_pct = TRUE) %>%
    filter(time <= tmax) %>%
    mutate(cuminc = 100 - median_pred,          # cumulative incidence (%)
           lo = 100 - ci_upper, hi = 100 - ci_lower,
           stratum = factor(stratum_labels[pred_set], levels = stratum_labels),
           label = relabel_arm(arm, trial1),
           is_control = grepl("^Control", label)) %>%
    filter(!is_control | label %in% keep_controls)

  int_labels <- ci %>% filter(!is_control) %>% pull(label) %>% unique()
  int_labels <- int_labels[order(as.numeric(sub("-.*$", "", int_labels)))]
  pal <- c(setNames(red_gradient(length(int_labels)), int_labels),
           "Control ≥ 24mo" = "grey35", "Control ≥ 60mo" = "grey55")
  ci$label <- factor(ci$label, levels = c(int_labels, keep_controls))

  ggplot(ci, aes(time, cuminc, group = label, colour = label)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = label), alpha = 0.12, colour = NA) +
    geom_line(aes(linetype = is_control), linewidth = 0.8) +
    facet_wrap(~ stratum, nrow = 1) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"), guide = "none") +
    scale_x_continuous(breaks = seq(0, tmax, 12)) +
    scale_y_continuous(limits = c(0, 40)) +
    labs(x = "Time since time zero (months)",
         y = "Cumulative incidence of breast cancer events (%)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

fig2 <- plot_cuminc(bootstrap_first, trial1 = TRUE,  keep_controls = "Control ≥ 24mo")
fig3 <- plot_cuminc(bootstrap_second, trial1 = FALSE)

ggsave(file.path(figures_dir, "Figure2_trial1_cuminc.pdf"), fig2, width = 30, height = 11, units = "cm")
ggsave(file.path(figures_dir, "Figure3_trial2_cuminc.pdf"), fig3, width = 30, height = 11, units = "cm")

################################################################################
# Forest plots of 5-year risk differences (Figure 4)
################################################################################

# One forest panel (points + 95% CI whiskers + value label) faceted by stratum,
# for a given control arm.
plot_forest <- function(df_diff, control_label, trial1 = FALSE) {

  d <- df_diff %>%
    filter(control_arm == control_label) %>%
    mutate(window = if (trial1) "18-30mo" else as.character(interruption_window)) %>%
    mutate(window = factor(window,
             levels = if (trial1) "18-30mo"
                      else c("6-12","12-18","18-24","24-30","30-36",
                             "36-42","42-48","48-54","54-60")),
           dataset = factor(dataset,
             levels = c("Overall", "Low/Intermediate risk", "High risk")))

  pal <- red_gradient(nlevels(d$window)); names(pal) <- levels(d$window)

  ggplot(d, aes(window, median_diff_pct, colour = window)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_linerange(aes(ymin = diff_lower_pct, ymax = diff_upper_pct), linewidth = 0.7) +
    geom_point(size = 2.6) +
    geom_text(aes(y = diff_upper_pct, label = sprintf("%.1f", median_diff_pct)),
              vjust = -0.6, size = 3, show.legend = FALSE) +
    facet_wrap(~ dataset, nrow = 1) +
    scale_colour_manual(values = pal, guide = "none") +
    labs(x = "Interruption window (months of ET)",
         y = "5-year difference (percentage points)",
         title = control_label) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 0, size = 8, face = "bold"),
          panel.border = element_rect(fill = NA, colour = "grey70"))
}

diff2 <- extract_5year_difference(bootstrap_second, target_time = 60)$diff_vs_control %>%
  mutate(interruption_window = sub("^x([0-9]+-[0-9]+)_.*$", "\\1", test_arm),
         control_arm = ifelse(grepl("24", control_arm), "Control ≥ 24mo", "Control ≥ 60mo"),
         dataset = factor(dataset, levels = c("overall","low_risk","high_risk"),
                          labels = c("Overall","Low/Intermediate risk","High risk")))

fig4_A <- plot_forest(diff2, "Control ≥ 24mo")   # Figure 4A
fig4_B <- plot_forest(diff2, "Control ≥ 60mo")   # Figure 4B
fig4 <- fig4_A / fig4_B                                # stack A over B (patchwork)

ggsave(file.path(figures_dir, "Figure4_trial2_forest.pdf"), fig4, width = 34, height = 20, units = "cm")

cat("\nFigures written to ", figures_dir, "\n")

################################################################################
# Pregnancy-seeking panels (Figures 2B, 3B, 4C-D)
################################################################################
# Re-run the pregnancy-seeking sensitivity analysis (censoring individuals
# without pregnancy-planning indicators; see Supplementary Methods), obtain its
# bootstrap object, and call plot_cuminc()/plot_forest() on it exactly as above.
