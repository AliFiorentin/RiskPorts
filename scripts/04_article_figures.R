## =============================================================================
## 04_article_figures.R -- Covariate balance figures and appendix table
##
## Produces, from data/baserisk1.dta:
##   figures/Figure_1  covariate balance, standardized mean differences and
##                     variance ratios, for the three matching algorithms
##   figures/Figure_2  common support, propensity score before and after
##                     nearest-neighbour matching
##   tables/table_a1_balance.tex  group means behind the balance analysis
##
## Each figure is written as vector PDF (used by LaTeX) and 600 dpi PNG.
##
## This script does NOT change the estimation. It reproduces the matching calls
## of 02_psm_algorithms.R exactly -- same formulas, same methods, same
## estimands, no caliper -- and changes only how results are presented:
## English labels, high resolution, vector output, colour-blind safe palette.
##
## Divergences between the code and the text of the manuscript are recorded in
## comments below, for the authors to decide on. They are not "fixed" here:
##
##   (a) The manuscript states "a caliper of 0.001 was applied". No matchit()
##       call in the project uses a caliper. Kept without one, as in the code.
##   (b) Full matching uses estimand = "ATE", while Equation (2) of the
##       manuscript defines the ATT. Kept as "ATE", as in the code.
##   (c) The third algorithm is labelled "radius" in the scripts and "Radius
##       Matching" in Table A.1, but it is method = "optimal", ratio = 2.
##
## Unlike 01 and 02, this script is self-contained: it can be run from start to
## finish with source(), from any working directory.
##
## Encoding: ASCII only, to avoid locale problems between Windows and Linux.
## =============================================================================

## -----------------------------------------------------------------------------
## 1. Configuration
## -----------------------------------------------------------------------------

## Repository root, discovered automatically. Works under Rscript, source() and
## the RStudio Source button, from any working directory. Override only if the
## data lives outside the repository.
find_root <- function() {
  marker <- file.path("data", "baserisk1.dta")

  climb <- function(from) {
    dir <- normalizePath(from, winslash = "/", mustWork = FALSE)
    for (i in 1:6) {
      if (file.exists(file.path(dir, marker))) return(dir)
      parent <- dirname(dir)
      if (identical(parent, dir)) break
      dir <- parent
    }
    NULL
  }

  ## (a) Rscript passes the path in --file=
  args <- commandArgs(trailingOnly = FALSE)
  arg_file <- grep("^--file=", args, value = TRUE)
  if (length(arg_file) == 1) {
    root <- climb(dirname(sub("^--file=", "", arg_file)))
    if (!is.null(root)) return(root)
  }

  ## (b) source() records the file path in the calling frame
  for (n in seq_len(sys.nframe())) {
    ofile <- tryCatch(sys.frame(n)$ofile, error = function(e) NULL)
    if (!is.null(ofile)) {
      root <- climb(dirname(ofile))
      if (!is.null(root)) return(root)
    }
  }

  ## (c) RStudio: path of the document open in the editor
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
    path <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                     error = function(e) NULL)
    if (!is.null(path) && nzchar(path)) {
      root <- climb(dirname(path))
      if (!is.null(root)) return(root)
    }
  }

  ## (d) Last resort: climb from the working directory
  root <- climb(getwd())
  if (!is.null(root)) return(root)

  stop("Could not locate the repository root (the folder containing ", marker, ").\n",
       "Working directory: ", getwd(), "\n",
       "Set the path manually, for example:\n",
       "  PROJECT_ROOT <- \"C:/path/to/RiskPorts\"", call. = FALSE)
}

PROJECT_ROOT <- find_root()

DATA_FILE  <- file.path(PROJECT_ROOT, "data", "baserisk1.dta")
FIGURE_DIR <- file.path(PROJECT_ROOT, "figures")
TABLE_DIR  <- file.path(PROJECT_ROOT, "tables")

message("Repository root: ", PROJECT_ROOT)

FIRST_YEAR <- 2018      # sample window used in the paper
SEED       <- 20260728  # matching with replacement breaks ties at random

## Output size, in inches (Elsevier full page width is about 7.5 in).
FIG_WIDTH  <- 7.5
FIG_HEIGHT <- 5.0
FIG_DPI    <- 600

## -----------------------------------------------------------------------------
## 2. Packages
## -----------------------------------------------------------------------------

packages <- c("haven", "MatchIt", "cobalt", "ggplot2")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nInstall with: install.packages(c(\"",
       paste(missing, collapse = "\", \""), "\"))", call. = FALSE)
}
## optmatch is required by method = "full" and method = "optimal"; checked
## separately because MatchIt's error in that case is not informative.
if (!requireNamespace("optmatch", quietly = TRUE)) {
  stop("Full and optimal matching require optmatch: install.packages(\"optmatch\")",
       call. = FALSE)
}

library(ggplot2)

## -----------------------------------------------------------------------------
## 3. Data
## -----------------------------------------------------------------------------

if (!file.exists(DATA_FILE)) {
  stop("Dataset not found: ", DATA_FILE, call. = FALSE)
}

base <- haven::read_dta(DATA_FILE)

covariates <- c("movimentacao", "terminais",
                "d_carga1", "d_carga2", "d_carga3", "d_carga4", "p_publico")
required <- c("parada", "ano", covariates)
absent <- setdiff(required, names(base))
if (length(absent) > 0) {
  stop("Columns missing from ", basename(DATA_FILE), ": ",
       paste(absent, collapse = ", "), call. = FALSE)
}

## Year filter (a no-op when the file is already restricted) and squared cargo
## volume. If the column already exists it is preserved, not recomputed.
sample_df <- base[base$ano >= FIRST_YEAR, ]
if (!"movimentacao2" %in% names(sample_df)) {
  sample_df$movimentacao2 <- sample_df$movimentacao^2
}

## The two formulas of 02_psm_algorithms.R, reproduced as they are. The original
## script uses `terminais` only in the third algorithm; this is not harmonised.
model_without_terminals <- parada ~ movimentacao + movimentacao2 +
  d_carga1 + d_carga2 + d_carga3 + d_carga4 + p_publico

model_with_terminals <- parada ~ movimentacao + movimentacao2 + terminais +
  d_carga1 + d_carga2 + d_carga3 + d_carga4 + p_publico

## All algorithms must run on EXACTLY the same rows, otherwise the weight
## vectors do not align in the cobalt comparison.
sample_df <- sample_df[stats::complete.cases(sample_df[, all.vars(model_with_terminals)]), ]
sample_df$parada <- as.numeric(sample_df$parada)

message("Analysis sample: ", nrow(sample_df), " vessel calls (",
        sum(sample_df$parada == 1), " treated / ",
        sum(sample_df$parada == 0), " control)")

## -----------------------------------------------------------------------------
## 4. Matching -- the three algorithms of Appendix A
## -----------------------------------------------------------------------------

set.seed(SEED)

## (i) Nearest neighbour with replacement -- 02_psm_algorithms.R, line 4.
m_nearest <- MatchIt::matchit(model_without_terminals, data = sample_df,
                              method = "nearest", m.order = "largest",
                              replace = TRUE)

## (ii) Optimal full matching -- 02_psm_algorithms.R, line 25.
m_optimal <- MatchIt::matchit(model_without_terminals, data = sample_df,
                              distance = "glm", method = "full",
                              estimand = "ATE")

## (iii) Third algorithm -- 02_psm_algorithms.R, line 34. The only one that
##       includes `terminais`, and the only one labelled "radius" despite being
##       method = "optimal", ratio = 2. Reproduced as it is.
m_radius <- MatchIt::matchit(model_with_terminals, data = sample_df,
                             distance = "glm", method = "optimal", ratio = 2)

## -----------------------------------------------------------------------------
## 5. Labels
## -----------------------------------------------------------------------------

## Cargo dummy mapping verified by cross-tabulation against `perfil_carga`
## (exact 1:1 correspondence). It does not follow the order in which the
## categories are listed in the paper. See data/README.md.
labels <- c(
  distance      = "Propensity score",
  movimentacao  = "Cargo volume (tons)",
  movimentacao2 = "Cargo volume, squared",
  terminais     = "Terminal",
  d_carga1      = "Containerised cargo",
  d_carga2      = "General cargo",
  d_carga3      = "Liquid and gas bulk",
  d_carga4      = "Dry bulk",
  p_publico     = "Public port"
)

## cobalt accepts a list of matchit objects in `weights` and extracts the
## weights from each, allowing the three algorithms to share one plot.
##
## The formula below only decides WHICH COVARIATES ARE DISPLAYED -- it is the
## union of the two, so that `terminais` appears. It changes no estimation: the
## weights come from the matchit objects, each fitted with its own formula.
balance <- cobalt::bal.tab(
  model_with_terminals, data = sample_df,
  weights   = list("Nearest neighbour" = m_nearest,
                   "Optimal"           = m_optimal,
                   "Radius"            = m_radius),
  stats     = c("mean.diffs", "variance.ratios"),
  s.d.denom = "treated", un = TRUE,
  distance  = data.frame(distance = m_nearest$distance)
)

article_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    plot.title       = element_text(face = "bold", size = 11)
  )

## Okabe-Ito palette, safe for colour vision deficiency, as required by the
## Marine Policy artwork guidelines. Shapes vary alongside colour so the plot
## remains readable in greyscale.
COLOURS <- c("#000000", "#E69F00", "#56B4E9", "#009E73")
SHAPES  <- c(16, 17, 15, 18)

## -----------------------------------------------------------------------------
## 6. Figures
## -----------------------------------------------------------------------------

## With two panels, love.plot returns a gtable rather than a ggplot. Adding
## "+ theme" to that object fails silently and yields a blank figure, so the
## theme is set globally before the call.
theme_set(article_theme)

## A single plot with both panels. Two reasons not to split them:
##   - with a shared y axis, ALL covariates stay listed in both panels; split
##     apart, the dummies vanish from the variance panel, because variance
##     ratios are not defined for binary variables;
##   - `var.order` stays at its default, preserving the formula order, which is
##     the order used in the paper.
## stars = "raw" marks with an asterisk the covariates displayed as raw
## differences in proportions, matching the original figure.
fig_balance <- cobalt::love.plot(
  balance, stats = c("mean.diffs", "variance.ratios"), abs = FALSE,
  var.names = labels, thresholds = c(m = 0.1, v = 2), stars = "raw",
  drop.distance = FALSE, title = NULL,
  colors = COLOURS, shapes = SHAPES
)

fig_support <- cobalt::bal.plot(
  m_nearest, var.name = "distance", which = "both", type = "density",
  sample.names = c("Before matching", "After matching")
) +
  labs(title = NULL, x = "Propensity score", y = "Density") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"),
                    labels = c("Control (no stoppage)", "Treated (stoppage)")) +
  article_theme

## -----------------------------------------------------------------------------
## 7. Export -- vector PDF plus 600 dpi PNG
## -----------------------------------------------------------------------------

if (!dir.exists(FIGURE_DIR)) dir.create(FIGURE_DIR, recursive = TRUE)

## Accepts both ggplot and gtable (love.plot with two panels returns a gtable).
save_figure <- function(plot, name, width = FIG_WIDTH, height = FIG_HEIGHT) {
  pdf_out <- file.path(FIGURE_DIR, paste0(name, ".pdf"))
  png_out <- file.path(FIGURE_DIR, paste0(name, ".png"))
  draw <- function(device, file, ...) {
    device(file, width = width, height = height, ...)
    on.exit(grDevices::dev.off())
    grid::grid.newpage(); grid::grid.draw(plot)
  }
  if (inherits(plot, "ggplot")) {
    ggsave(pdf_out, plot, width = width, height = height, units = "in")
    ggsave(png_out, plot, width = width, height = height, units = "in", dpi = FIG_DPI)
  } else {
    draw(grDevices::pdf, pdf_out)
    draw(grDevices::png, png_out, units = "in", res = FIG_DPI)
  }
  message("written: ", pdf_out)
  message("written: ", png_out)
  invisible(NULL)
}

save_figure(fig_balance, "Figure_1", width = 9.0, height = 5.5)
save_figure(fig_support, "Figure_2", height = 4.0)

## -----------------------------------------------------------------------------
## 8. Table A.1 -- built from the SAME matchit objects as the figures
##
## The appendix balance table is generated here rather than transcribed by
## hand, so that table and figure can never diverge. The manuscript reads the
## .tex file through \input{}.
## -----------------------------------------------------------------------------

if (!dir.exists(TABLE_DIR)) dir.create(TABLE_DIR, recursive = TRUE)

## Group means, before and after matching, for one matchit object.
balance_means <- function(m) {
  b <- cobalt::bal.tab(m, un = TRUE, disp = "means", stats = "mean.diffs",
                       s.d.denom = "treated")
  columns <- c("M.1.Un", "M.0.Un", "M.1.Adj", "M.0.Adj")
  absent <- setdiff(columns, names(b$Balance))
  if (length(absent) > 0) {
    stop("cobalt::bal.tab did not return the columns ", paste(absent, collapse = ", "),
         ".\ncobalt version: ", as.character(utils::packageVersion("cobalt")),
         call. = FALSE)
  }
  list(means = b$Balance[, columns, drop = FALSE], n = b$Observations)
}

## Thousands separator above 1,000; four decimals below.
fmt <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) return("--")
    if (abs(v) >= 1000) formatC(v, format = "f", digits = 2, big.mark = ",")
    else formatC(v, format = "f", digits = 4)
  }, character(1))
}

algorithms <- list("Nearest-Neighbor Matching" = m_nearest,
                   "Optimal Matching"          = m_optimal,
                   "Radius Matching"           = m_radius)

balances <- lapply(algorithms, balance_means)

## Table body: one row per covariate, four columns per algorithm.
variables <- rownames(balances[[1]]$means)
body <- vapply(variables, function(v) {
  cells <- unlist(lapply(balances, function(b) fmt(unlist(b$means[v, ]))))
  label <- if (v %in% names(labels)) labels[[v]] else v
  paste0(label, " & ", paste(cells, collapse = " & "), " \\\\")
}, character(1))

## Header sample sizes, read from the objects themselves.
header_n <- unlist(lapply(balances, function(b) {
  n <- b$n
  c(sprintf("\\makecell{Unmatched\\\\Treated\\\\(n=%s)", format(round(n["All", "Treated"]), big.mark = ",")),
    sprintf("\\makecell{Unmatched\\\\Control\\\\(n=%s)", format(round(n["All", "Control"]), big.mark = ",")),
    sprintf("\\makecell{Matched\\\\Treated\\\\(n=%s)",   format(round(n[2, "Treated"]),     big.mark = ",")),
    sprintf("\\makecell{Matched\\\\Control\\\\(n=%s)",   format(round(n[2, "Control"]),     big.mark = ",")))
}))
header_n <- paste0(header_n, "}")

n_alg <- length(algorithms)
groups <- paste(sprintf("\\multicolumn{4}{c}{\\textbf{%s}}", names(algorithms)), collapse = " & ")
cmidrules <- paste(vapply(seq_len(n_alg), function(i)
  sprintf("\\cmidrule(lr){%d-%d}", 2 + (i - 1) * 4, 1 + i * 4), character(1)), collapse = "")

table_tex <- c(
  "%% GENERATED BY scripts/04_article_figures.R -- DO NOT EDIT BY HAND.",
  "%% Re-run the script to update.",
  paste0("\\begin{tabular}{l", strrep("c", n_alg * 4), "}"),
  "\\toprule",
  paste0(" & ", groups, " \\\\"),
  cmidrules,
  paste0("\\textbf{Variable} & ", paste(header_n, collapse = " & "), " \\\\"),
  "\\midrule",
  unname(body),
  "\\bottomrule",
  "\\end{tabular}"
)

tex_file <- file.path(TABLE_DIR, "table_a1_balance.tex")
writeLines(table_tex, tex_file)
message("written: ", tex_file)

## Same numbers as CSV, for checking outside LaTeX.
csv <- do.call(rbind, lapply(names(balances), function(name) {
  m <- balances[[name]]$means
  data.frame(algorithm = name, variable = rownames(m), m, row.names = NULL,
             check.names = FALSE)
}))
csv_file <- file.path(TABLE_DIR, "table_a1_balance.csv")
utils::write.csv(csv, csv_file, row.names = FALSE)
message("written: ", csv_file)

message("\nDone. Figures are referenced without an extension in the manuscript, ",
        "so pdfLaTeX picks up the .pdf files automatically.")
