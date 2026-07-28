# RiskPorts

Data and code for **"Climate Variables and the Costs Generated for the Port Sector: The Case of a Port in Southern Brazil"**, submitted to *Marine Policy*.

The study estimates the operational delays and the opportunity cost caused by weather-related stoppages at the São Francisco do Sul port complex (Santa Catarina, Brazil), combining port operations records with meteorological and wave data. Identification relies on **Propensity Score Matching** followed by regression with meteorological controls.

## What is here

```
data/
  baserisk1.dta        analysis dataset: 4,332 vessel calls, 2018-2021
  README.md            codebook and provenance
scripts/
  01_matching_psm.R    propensity score, matching, and all regressions
  02_psm_algorithms.R  covariate balance for the three matching algorithms
  03_final_analysis.do logit/probit models and marginal effects (Stata)
  04_article_figures.R publication figures and the appendix balance table
```

## Reproducing the results

The analysis runs in **R** (tested on 4.6.1); one script requires **Stata**.

```r
install.packages(c("haven", "MatchIt", "cobalt", "optmatch", "Matching",
                   "rbounds", "ggplot2", "lmtest", "sandwich", "dplyr"))
```

Run in order. Scripts `01` and `02` reproduce the estimates and the balance
diagnostics; `04` regenerates the figures and the appendix table as they appear
in the paper.

```r
source("scripts/01_matching_psm.R")
source("scripts/02_psm_algorithms.R")
source("scripts/04_article_figures.R")
```

`04_article_figures.R` is self-contained and can be run on its own.

### One caveat on exact reproduction

Nearest-neighbour matching **with replacement** breaks ties at random, and the
original analysis did not set a seed. Point estimates in Table 2 of the paper
reproduce exactly; the Rosenbaum bounds in Table 3 are sensitive to the draw
(at $\Gamma = 1.6$ the bound ranges from 0.71 to 0.99 across seeds). Script `04`
sets a seed so that its output is stable.

## Data sources and availability

**Port operations** — originally extracted from the public waterway statistics
of the Brazilian National Waterway Transportation Agency (ANTAQ),
<https://www.gov.br/antaq/>. **This extraction is no longer reproducible from
the source:** ANTAQ has restructured its statistical system and discontinued
publication of the port revenue variable used in the opportunity cost
estimates. The dataset deposited here preserves the version used in the
analysis.

**Meteorological series** — provided by the Center for Information on Water and
Meteorological Resources of Santa Catarina (EPAGRI/CIRAM),
<https://ciram.epagri.sc.gov.br/>. The values matched to each vessel call are
included in `data/baserisk1.dta`. The original station-level series are not
redistributed here and may be obtained from EPAGRI/CIRAM upon request.

**Wave data** — significant wave height simulated under the "Riskports" project.

## Citation

Please cite the article. If you use the dataset or the code directly, please
also cite this repository through its archived DOI.

## License

Code is released under the MIT License. The dataset is released under
CC BY 4.0. See `LICENSE`.
