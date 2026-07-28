## =============================================================================
## 02_psm_algorithms.R -- Covariate balance across matching algorithms
##
## Compares the three propensity score matching algorithms reported in
## Appendix A of the paper: nearest neighbour, optimal (full) and the one
## labelled "radius".
##
## Originally written as a continuation of 01_matching_psm.R, relying on the
## objects it leaves in memory. The block below loads the data if they are not
## present, so the script can also be run on its own. No estimation was changed.
##
## Package names were corrected: the original called install.packages("opmatch")
## and install.packages("ggplot"), neither of which exists on CRAN -- the
## correct names are optmatch and ggplot2, so those lines always failed.
##
## Encoding: ASCII only.
## =============================================================================

library(MatchIt)
library(cobalt)
library(ggplot2)

## Continue from 01_matching_psm.R when its objects are in memory; otherwise
## rebuild the analysis sample from the dataset.
if (!exists("baserisk1")) {
  library(haven)
  root <- if (file.exists(file.path("data", "baserisk1.dta"))) "." else ".."
  baserisk1 <- read_dta(file.path(root, "data", "baserisk1.dta"))
  baserisk1 <- baserisk1[baserisk1$ano >= 2018, ]
  baserisk1$parada <- as.numeric(baserisk1$parada)
  if (!"movimentacao2" %in% names(baserisk1)) {
    baserisk1$movimentacao2 <- baserisk1$movimentacao^2
  }
  message("Loaded ", nrow(baserisk1), " vessel calls from data/baserisk1.dta")
}

## -----------------------------------------------------------------------------
## Nearest neighbour matching
## -----------------------------------------------------------------------------

mod_match <- matchit(parada ~ movimentacao + movimentacao2 + d_carga1 + d_carga2 +
                       d_carga3 + d_carga4 + p_publico,
                     data = baserisk1, method = "nearest", m.order = "largest",
                     replace = TRUE)
summary(mod_match)

## Covariate balance
love.plot(mod_match, binary = "raw",
          stars = "raw",
          stat = c("m", "v"),
          grid = TRUE,
          thresholds = c(m = 5, v = 10))

hist(mod_match$distance, xlab = "Distance", xlim = c(0, 1.5),
     main = "Difference between groups", probability = TRUE)
lines(density(mod_match$distance))

## -----------------------------------------------------------------------------
## Optimal (full) matching
## -----------------------------------------------------------------------------

mod_match1 <- matchit(parada ~ movimentacao + movimentacao2 + d_carga1 + d_carga2 +
                        d_carga3 + d_carga4 + p_publico,
                      data = baserisk1, distance = "glm", method = "full",
                      estimand = "ATE")
## The original called summary(mod_match) here, which summarises the previous
## object. Corrected to the object just fitted.
summary(mod_match1)

love.plot(mod_match1, binary = "raw",
          stars = "raw",
          stat = c("m", "v"),
          grid = TRUE,
          thresholds = c(m = 5, v = 10))

## -----------------------------------------------------------------------------
## "Radius" matching
##
## NOTE: the call is method = "optimal" with ratio = 2, that is, optimal 1:2
## pair matching -- not radius (caliper) matching. The label is kept as in the
## original because this is the call that produced Table A.1 of the paper.
## It is also the only one of the three that includes `terminais`.
## -----------------------------------------------------------------------------

mod_match2 <- matchit(parada ~ movimentacao + movimentacao2 + terminais +
                        d_carga1 + d_carga2 + d_carga3 + d_carga4 + p_publico,
                      data = baserisk1, distance = "glm", method = "optimal",
                      ratio = 2)
summary(mod_match2)

love.plot(mod_match2, binary = "raw",
          stars = "raw",
          stat = c("m", "v"),
          grid = TRUE,
          thresholds = c(m = 5, v = 10))

## Note: this script draws to the screen only. The figures published with the
## paper are produced by 04_article_figures.R, which writes them to disk.
