################################################################################
# 00_setup.R
# Paths, packages and plotting theme shared by all analysis scripts.
#
# Target trial emulation of endocrine therapy (ET) interruption for pregnancy
# in HR+ early breast cancer (SNDS). See README.md.
#
# EDIT the paths below to point to your local copy of the (restricted) SNDS-derived
# datasets. The SNDS data are not public (CNIL authorisation 922143v1); see README.
################################################################################

## ---- Paths -----------------------------------------------------------------
# Directory holding the person-month .RData datasets (not shared publicly).
data_dir    <- Sys.getenv("POSITIVE_DATA_DIR", unset = "data/")
results_dir <- "results/"
figures_dir <- "figures/"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Packages --------------------------------------------------------------

pkgs <- c("data.table", "dplyr", "tidyr", "purrr",
          "splines", "speedglm", "broom", "ggplot2", "patchwork")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nInstall them with install.packages().")
}
invisible(lapply(pkgs, library, character.only = TRUE))

## ---- Plotting theme --------------------------------------------------------
ggplot2::theme_set(ggplot2::theme_minimal())

## Reproducibility: bootstrap scripts set their own seed with set.seed().
message("Setup complete. data_dir = ", normalizePath(data_dir, mustWork = FALSE))
