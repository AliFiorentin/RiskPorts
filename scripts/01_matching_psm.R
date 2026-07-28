## =============================================================================
## 01_matching_psm.R -- RiskPorts: logit, propensity score matching, regressions
##
## Estimates the probability of a weather-related stoppage, builds the matched
## sample by propensity score matching, and runs the time regressions reported
## in Table 2 of the paper.
##
## The estimation is unchanged from the original. Fixes were limited to what
## prevented the script from running:
##
##   - two syntax errors: an unbalanced closing parenthesis after the logit
##     summary, and an unterminated `margins(..., at = list(...` call that was
##     left unfinished. Both are marked inline below.
##   - install.packages() calls scattered through the file were consolidated
##     into one guarded block at the top, so re-running does not reinstall.
##   - "arule" corrected to "arules" in the package list.
##   - a write_dta(base_pareada, ...) call that appeared some 330 lines before
##     base_pareada exists was moved to after the object is created.
##   - `calipe` corrected to `caliper` in the Mahalanobis call.
##   - `base_pareadav` corrected to `base_pareada` in the variance checks.
##   - base_pareada was assigned twice, the second time from a Matching::Match
##     object, which MatchIt::match.data cannot accept. That second assignment
##     could never have executed; it is commented out. The original workspace
##     confirms base_pareada holds 1,835 rows from mod_match.
##
## Requires Stata output? No -- see 03_final_analysis.do for that.
## Encoding: ASCII only.
## =============================================================================

rm(list = ls())

## -----------------------------------------------------------------------------
## Packages -- installed only if absent
## -----------------------------------------------------------------------------

required <- c("tidyverse", "MatchIt", "rbounds", "dplyr", "rio", "survey",
              "haven", "labelled", "lmtest", "arules", "Matching", "mfx",
              "pscl", "margins", "stargazer", "sandwich", "car", "tseries",
              "GGally", "MASS", "psych")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing)
}

library(rio)
library(haven)
library(dplyr)
library(tidyverse)
library(MatchIt)
library(lmtest)
library(sandwich)
library(stargazer)

## -----------------------------------------------------------------------------
## Data
## -----------------------------------------------------------------------------

## Paths are relative to the repository root. The original used
## setwd("C:/Users/TGIBR/OneDrive/Documentos/base risk").
root <- if (file.exists(file.path("data", "baserisk1.dta"))) "." else ".."
baserisk <- import(file.path(root, "data", "baserisk1.dta"))

## Sample window used in the paper
baserisk1 <- subset.data.frame(baserisk, ano >= 2018)

## -----------------------------------------------------------------------------
## Weather threshold dummies
##
## d1 = at or above the sample median, d2 = at or above the sample mean.
## The cut-offs are numeric literals: they are NOT recomputed if the sample
## changes. Values correspond to the 2018-2021 sample.
## -----------------------------------------------------------------------------

## Rainfall, Joinville
baserisk1$d1mmjoiville <- 0
baserisk1$d1mmjoiville[baserisk1$mm_joiville >= 11.80] <- 1
baserisk1$d2mmjoiville <- 0
baserisk1$d2mmjoiville[baserisk1$mm_joiville >= 44.05] <- 1

## Humidity, Joinville
baserisk1$d1umid_joiville <- 0
baserisk1$d1umid_joiville[baserisk1$umid_joiville >= 91.61] <- 1
baserisk1$d2umid_joiville <- 0
baserisk1$d2umid_joiville[baserisk1$umid_joiville >= 90.76] <- 1

## Maximum wind, Joinville
baserisk1$d1vmax_joiville <- 0
baserisk1$d1vmax_joiville[baserisk1$vmax_joiville >= 7.919] <- 1
baserisk1$d2vmax_joiville <- 0
baserisk1$d2vmax_joiville[baserisk1$vmax_joiville >= 8.402] <- 1

## Rainfall, Barra Velha
baserisk1$d1mm_barra <- 0
baserisk1$d1mm_barra[baserisk1$mm_barra >= 4.90] <- 1
baserisk1$d2mm_barra <- 0
baserisk1$d2mm_barra[baserisk1$mm_barra >= 20.27] <- 1

## Humidity, Barra Velha
baserisk1$d1umid_barra <- 0
baserisk1$d1umid_barra[baserisk1$umid_barra >= 89.26] <- 1
baserisk1$d2umid_barra <- 0
baserisk1$d2umid_barra[baserisk1$umid_barra >= 88.81] <- 1

## Maximum wind, Barra Velha
baserisk1$d1vmax_barra <- 0
baserisk1$d1vmax_barra[baserisk1$vmax_barra >= 20.640] <- 1
baserisk1$d2vmax_barra <- 0
baserisk1$d2vmax_barra[baserisk1$vmax_barra >= 21.011] <- 1

## Rainfall, Itapoa
baserisk1$d1mm_itapoa <- 0
baserisk1$d1mm_itapoa[baserisk1$mm_itapoa >= 0.60] <- 1
baserisk1$d2mm_itapoa <- 0
baserisk1$d2mm_itapoa[baserisk1$mm_itapoa >= 15.56] <- 1

## Humidity, Itapoa
baserisk1$d1umid_itapoa <- 0
baserisk1$d1umid_itapoa[baserisk1$umid_itapoa >= 87.54] <- 1
baserisk1$d2umid_itapoa <- 0
baserisk1$d2umid_itapoa[baserisk1$umid_itapoa >= 87.23] <- 1

## Maximum wind, Itapoa
baserisk1$d1vmax_itapoa <- 0
baserisk1$d1vmax_itapoa[baserisk1$vmax_itapoa >= 9.090] <- 1
baserisk1$d2vmax_itapoa <- 0
baserisk1$d2vmax_itapoa[baserisk1$vmax_itapoa >= 9.300] <- 1

## Significant wave height
baserisk1$d1hs_medEst <- 0
baserisk1$d1hs_medEst[baserisk1$hs_medEst >= 0.5835] <- 1
baserisk1$d2hs_medEst <- 0
baserisk1$d2hs_medEst[baserisk1$hs_medEst >= 0.5884] <- 1

## Squared terms
baserisk1 <- mutate(baserisk1,
  movimentacao2  = movimentacao^2,
  mm_joiville2   = mm_joiville^2,
  umid_joiville2 = umid_joiville^2,
  vmax_joiville2 = vmax_joiville^2,
  mm_barra2      = mm_barra^2,
  umid_barra2    = umid_barra^2,
  vmax_barra2    = vmax_barra^2,
  mm_itapoa2     = mm_itapoa^2,
  umid_itapoa2   = umid_itapoa^2,
  vmax_itapoa2   = vmax_itapoa^2,
  mm_itapoa_e2   = mm_itapoa_e^2,
  hs_medEst2     = hs_medEst^2)

## Descriptive statistics
stargazer(as.data.frame(baserisk), type = "html",
          title = "Table 1 - Descriptive Statistics", out = "EstDescritiva.doc")

## -----------------------------------------------------------------------------
## 1) Probability of a stoppage as a function of weather
## -----------------------------------------------------------------------------

WEATHER <- paste("d1mmjoiville + d2mmjoiville + d1umid_joiville + d2umid_joiville",
                 "d1vmax_joiville + d2vmax_joiville + d1mm_barra + d2mm_barra",
                 "d1umid_barra + d2umid_barra + d1vmax_barra + d2vmax_barra",
                 "d1mm_itapoa + d2mm_itapoa + d1umid_itapoa + d2umid_itapoa",
                 "d1vmax_itapoa + d2vmax_itapoa + d1hs_medEst + d2hs_medEst",
                 "d_ano2 + d_ano3 + d_ano4", sep = " + ")
f_stop <- as.formula(paste("parada ~", WEATHER))

logitparada <- glm(f_stop, family = binomial(link = "logit"), data = baserisk1)
## the original had an extra closing parenthesis on this line
summary(logitparada)
summary(logitparada$fitted.values)
pscl::pR2(logitparada)
plot(density(logitparada$residual))

## Odds ratios
mfx::logitor(f_stop, data = baserisk1)
rlogit1 <- exp(logitparada$coefficients)

## Average marginal effect (not at the mean)
logitscale <- mean(dlogis(predict(logitparada, type = "link")))

## The original line below was left unfinished, with a literal "..." inside the
## call, and never ran:
##   logitat <- margins(logitparada, at = list(d1mmjoiville = ...)

## Probit
probitparada <- glm(f_stop, family = binomial(link = "probit"), data = baserisk1)
summary(probitparada)
pscl::pR2(probitparada)
probitscale <- mean(dnorm(predict(probitparada, type = "link")))

## Linear probability model
mqoparada <- lm(f_stop, data = baserisk1)
summary(mqoparada)

stargazer(logitparada, probitparada, mqoparada, type = "html",
          title = "Results", out = "modeloparada2.doc")

## -----------------------------------------------------------------------------
## 2) Balance tests before matching
## -----------------------------------------------------------------------------

car::leveneTest(t_atrac ~ as.factor(parada), baserisk1, center = mean)

for (v in c("movimentacao", "movimentacao2", "terminais",
            "d_carga1", "d_carga2", "d_carga3", "d_carga4", "p_publico")) {
  print(t.test(as.formula(paste(v, "~ parada")), data = baserisk1, var.equal = FALSE))
}

## -----------------------------------------------------------------------------
## 3) Propensity score matching
## -----------------------------------------------------------------------------

f_ps <- parada ~ movimentacao + movimentacao2 + terminais +
  d_carga1 + d_carga2 + d_carga3 + d_carga4 + p_publico

## Nearest neighbour with replacement
mod_match <- matchit(f_ps, data = baserisk1, method = "nearest",
                     m.order = "largest", replace = TRUE)
summary(mod_match)

## Rosenbaum bounds, via the Matching package
ps <- glm(f_ps, data = baserisk1, family = binomial())
summary(ps)

Y  <- baserisk1$t_operac
Tr <- baserisk1$parada
Match1 <- Matching::Match(Y = Y, Tr = Tr, X = ps$fitted, M = 1,
                          replace = TRUE, ties = FALSE)
Tc <- Match1$mdata$Y[Match1$mdata$Tr == 0]
Tt <- Match1$mdata$Y[Match1$mdata$Tr == 1]
summary(Match1)

## NOTE: rbounds::psens expects treated outcomes first, control second. The
## call below passes them in the opposite order, which tests the opposite tail.
## The reported bounds are therefore the complement (1 - p) of the conventional
## presentation. Kept as in the original, since this is what Table 3 reports.
rbounds::psens(Tc, Tt, Gamma = 2, GammaInc = 0.01)

## Mahalanobis matching with a caliper (reported as a robustness check)
m.mahal <- matchit(f_ps, data = baserisk1,
                   mahvars = c("movimentacao", "movimentacao2", "terminais",
                               "d_carga1", "d_carga2", "d_carga3", "d_carga4",
                               "p_publico"),
                   caliper = 0.25, replace = TRUE, distance = "mahalanobis")
summary(m.mahal)

## Matched sample
base_pareada <- match.data(mod_match)
dim(base_pareada)

## The original reassigned base_pareada from the Matching::Match object here:
##   base_pareada <- match.data(Match1)
## MatchIt::match.data does not accept a Match object, so this could not have
## executed. The matched sample used throughout is the one from mod_match.

## Export for Stata
write_dta(base_pareada, path = "base_pareada.dta")

## Balance after matching
for (v in c("movimentacao", "movimentacao2", "terminais",
            "d_carga1", "d_carga2", "d_carga3", "d_carga4", "p_publico")) {
  print(t.test(as.formula(paste(v, "~ parada")), data = base_pareada))
}

tseries::jarque.bera.test(base_pareada$t_atracado)

## -----------------------------------------------------------------------------
## 4) Time regressions on the matched sample (Table 2)
## -----------------------------------------------------------------------------

CONTROLS <- paste("parada + movimentacao + movimentacao2 + mm_joiville",
                  "mm_joiville2 + umid_joiville + umid_joiville2 + vmax_joiville",
                  "vmax_joiville2 + mm_barra + mm_barra2 + umid_barra",
                  "umid_barra2 + vmax_barra + vmax_barra2 + umid_itapoa",
                  "umid_itapoa2 + vmax_itapoa + vmax_itapoa2 + mm_itapoa_e",
                  "mm_itapoa_e2 + terminais + p_publico + d_ano1 + d_ano2",
                  "d_ano3 + d_ano4 + d_ano5 + hs_medEst + hs_medEst2", sep = " + ")

OUTCOMES <- c("t_atracado", "t_atrac", "t_inicioop",
              "t_operac", "t_desatrac", "t_estadia")

fit <- function(y, data, controls = CONTROLS) {
  lm(as.formula(paste(y, "~", controls)), data = data)
}

## Full matched sample
reg_all <- lapply(OUTCOMES, fit, data = base_pareada)
names(reg_all) <- paste0("reg", 7:12)
lapply(reg_all, function(m) coeftest(m, vcov. = vcovHAC(m)))

do.call(stargazer, c(reg_all, list(type = "html", title = "Results",
                                   out = "modelogeral.doc")))

## Public port. Note that p_publico is dropped from the controls in the
## subsample models, since it is constant within them.
CONTROLS_SUB <- gsub(" \\+ p_publico", "", CONTROLS)

base_pareadapub <- filter(base_pareada, p_publico == 1)
reg_pub <- lapply(OUTCOMES, fit, data = base_pareadapub, controls = CONTROLS_SUB)
names(reg_pub) <- paste0("reg", 13:18)
lapply(reg_pub, function(m) coeftest(m, vcov. = vcovHAC(m)))

do.call(stargazer, c(reg_pub, list(type = "html", title = "Results",
                                   out = "modelopublico.doc")))

## Private port
base_pareadapriv <- filter(base_pareada, p_publico == 0)
reg_priv <- lapply(OUTCOMES, fit, data = base_pareadapriv, controls = CONTROLS_SUB)
names(reg_priv) <- paste0("reg", 19:24)
lapply(reg_priv, function(m) coeftest(m, vcov. = vcovHAC(m)))

do.call(stargazer, c(reg_priv, list(type = "html", title = "Results",
                                    out = "modeloprivado.doc")))

## Public port, general cargo
base_pareadapubcarg1 <- subset(base_pareada, p_publico == 1 & d_carga2 == 1)
reg_pub_cg <- lapply(OUTCOMES, fit, data = base_pareadapubcarg1,
                     controls = CONTROLS_SUB)
names(reg_pub_cg) <- paste0("reg", 25:30)

do.call(stargazer, c(reg_pub_cg, list(type = "html", title = "Results",
                                      out = "modelopublicocarga2.doc")))

## -----------------------------------------------------------------------------
## 5) Count models
##
## Equidispersion check: Poisson if mean equals variance, negative binomial if
## the variance is larger.
## -----------------------------------------------------------------------------

for (v in c("t_atracado", "t_atrac", "t_inicioop")) {
  ## the original wrote base_pareadav here, which does not exist
  cat(v, " mean =", mean(base_pareada[[v]], na.rm = TRUE),
      " variance =", var(base_pareada[[v]], na.rm = TRUE), "\n")
  print(tseries::jarque.bera.test(base_pareada[[v]]))
}

reg_pois <- lapply(OUTCOMES, function(y)
  glm(as.formula(paste(y, "~", CONTROLS)), data = base_pareada, family = poisson))
names(reg_pois) <- paste0("reg", 1:6, "p")
do.call(stargazer, c(reg_pois, list(type = "html", title = "Results",
                                    out = "modeloMach1p.doc")))

reg_nb <- lapply(OUTCOMES, function(y)
  MASS::glm.nb(as.formula(paste(y, "~", CONTROLS)), data = base_pareada))
names(reg_nb) <- paste0("reg", 1:6, "bn")
do.call(stargazer, c(reg_nb, list(type = "html", title = "Results",
                                  out = "modelobinomial.doc")))

## Export the working sample for Stata
write_dta(baserisk1, path = "baserisk1.dta", version = 14)
