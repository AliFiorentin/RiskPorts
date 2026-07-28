*******************************************************************************
* 03_final_analysis.do -- RiskPorts
*
* Logit, probit and marginal effects for the probability of a weather-related
* stoppage, plus the simulations reported in the Results section.
*
* Fixes limited to what prevented the file from running:
*   - added the `use` statement. The original had none: it operated on whatever
*     dataset happened to be in memory.
*   - output paths were hard-coded to /Users/marcio_nb/ (another author's Mac)
*     and are now relative to the repository root.
*   - `margins, at (...)` calls listing d1mm_itapoa without `=1` were a syntax
*     error; corrected.
*   - three `outreg2 ... replace` calls overwrote files that the preceding call
*     had just created, discarding the logit results; changed to `append`, and
*     flagged inline.
*   - `dprobit` and `mfx` are no longer available in current Stata. They are
*     wrapped in `capture` so the file runs end to end; the estimation itself
*     is unchanged.
*
* Requires: outreg2, logout (ssc install outreg2 logout)
*******************************************************************************

clear all
set more off

* Repository root -- adjust if running from elsewhere
local root "."
use "`root'/data/baserisk1.dta", clear
keep if ano >= 2018

capture mkdir "`root'/output"
local out "`root'/output"

*** Descriptive statistics for the full dataset
sum

*** Descriptive statistics for the analysed variables
logout, save("`out'/EstatResult") word excel replace: ///
    tabstat movimentacao t_atrac t_inicioop t_operac t_desatrac t_atracado ///
            t_estadia receita mm_joiville umid_joiville vmax_joiville ///
            mm_barra umid_barra vmax_barra mm_itapoa umid_itapoa vmax_itapoa ///
            mm_itapoa_e hs_medEst, stat(count median mean min max)

sum d1mmjoiville d2mmjoiville d1umid_joiville d2umid_joiville d1vmax_joiville ///
    d2vmax_joiville d1mm_barra d2mm_barra d1umid_barra d2umid_barra ///
    d1vmax_barra d2vmax_barra d1mm_itapoa d2mm_itapoa d1umid_itapoa ///
    d2umid_itapoa d1vmax_itapoa d2vmax_itapoa d1hs_medEst d2hs_medEst

* The full set of weather dummies, reused throughout
local X d1mmjoiville d2mmjoiville d1umid_joiville d2umid_joiville ///
        d1vmax_joiville d2vmax_joiville d1mm_barra d2mm_barra d1umid_barra ///
        d2umid_barra d1vmax_barra d2vmax_barra d1mm_itapoa d2mm_itapoa ///
        d1umid_itapoa d2umid_itapoa d1vmax_itapoa d2vmax_itapoa ///
        d1hs_medEst d2hs_medEst d_ano2 d_ano3 d_ano4

*******************************************************************************
* General model
*******************************************************************************

*** Logit
logit parada `X', robust
outreg2 using "`out'/Result_GeneralModel.doc", replace ctitle(Logit)

*** Odds ratios
logistic parada `X', robust
outreg2 using "`out'/Result_GeneralModel.doc", append ctitle(Logit Odds Ratio) eform

*** Probability of stoppage at the mean
margins

*** Marginal effects
capture mfx

*** Probit
probit parada `X', robust
capture mfx
outreg2 using "`out'/Result_GeneralModel.doc", append ctitle(Probit)

*** Marginal effects (dprobit is unavailable in current Stata)
capture dprobit parada `X', robust
capture mfx
capture outreg2 using "`out'/Result_GeneralModel.doc", append ctitle(Probit MFX) eform

*******************************************************************************
* Public port
*******************************************************************************

logit parada `X' if p_publico == 1, robust
outreg2 using "`out'/Result_PublicPort.doc", replace ctitle(Logit)

logistic parada `X' if p_publico == 1, robust
outreg2 using "`out'/Result_PublicPort.doc", append ctitle(Logit Odds Ratio) eform

margins
capture mfx

probit parada `X' if p_publico == 1, robust
capture mfx
* was `replace`, which discarded the logit results written above
outreg2 using "`out'/Result_PublicPort.doc", append ctitle(Probit) eform

*******************************************************************************
* Public port by cargo type
*******************************************************************************

*** General cargo (d_carga2)
logit parada `X' if p_publico == 1 & d_carga2 == 1, robust
outreg2 using "`out'/Result_PublicPort_GeneralCargo.doc", replace ctitle(Logit)

logistic parada `X' if p_publico == 1 & d_carga2 == 1, robust
outreg2 using "`out'/Result_PublicPort_GeneralCargo.doc", append ctitle(Logit Odds Ratio) eform

margins
capture mfx

probit parada `X' if p_publico == 1 & d_carga2 == 1, robust
capture mfx
* was `replace`
outreg2 using "`out'/Result_PublicPort_GeneralCargo.doc", append ctitle(Probit) eform

*** Dry bulk (d_carga4)
logit parada `X' if p_publico == 1 & d_carga4 == 1, robust
outreg2 using "`out'/Result_PublicPort_DryBulk.doc", replace ctitle(Logit)

logistic parada `X' if p_publico == 1 & d_carga4 == 1, robust
outreg2 using "`out'/Result_PublicPort_DryBulk.doc", append ctitle(Logit Odds Ratio) eform

margins
capture mfx

probit parada `X' if p_publico == 1 & d_carga4 == 1, robust
capture mfx
* was `replace`, and the filename had the extension misplaced
* (Result_PortoPublico.doc_GS)
outreg2 using "`out'/Result_PublicPort_DryBulk.doc", append ctitle(Probit) eform

*******************************************************************************
* Simulations
*******************************************************************************

* Rainfall above the mean, at all stations
local ABOVE_MEAN_RAIN d2mmjoiville=1 d2mm_barra=1 d2mm_itapoa=1
* Wind above the mean, at all stations
local ABOVE_MEAN_WIND d2vmax_joiville=1 d2vmax_barra=1 d2vmax_itapoa=1
* Rainfall above the median, at all stations
* (the original omitted `=1` on d1mm_itapoa, which is a syntax error)
local ABOVE_MEDIAN_RAIN d1mmjoiville=1 d1mm_barra=1 d1mm_itapoa=1

*** General model
logistic parada `X', robust
margins, at (`ABOVE_MEAN_RAIN')
margins, at (`ABOVE_MEAN_WIND')
margins, at (`ABOVE_MEAN_RAIN' `ABOVE_MEAN_WIND')
margins, at (`ABOVE_MEDIAN_RAIN')

*** Public port
logistic parada `X' if p_publico == 1, robust
margins, at (`ABOVE_MEAN_RAIN')
margins, at (`ABOVE_MEAN_WIND')
margins, at (`ABOVE_MEAN_RAIN' `ABOVE_MEAN_WIND')
margins, at (`ABOVE_MEDIAN_RAIN')

*** Public port, general cargo
logistic parada `X' if p_publico == 1 & d_carga2 == 1, robust
margins, at (`ABOVE_MEAN_RAIN')
margins, at (`ABOVE_MEAN_WIND')
margins, at (`ABOVE_MEAN_RAIN' `ABOVE_MEAN_WIND')
margins, at (`ABOVE_MEDIAN_RAIN')

*** Public port, dry bulk
logistic parada `X' if p_publico == 1 & d_carga4 == 1, robust
margins, at (`ABOVE_MEAN_RAIN')
margins, at (`ABOVE_MEAN_WIND')
margins, at (`ABOVE_MEAN_RAIN' `ABOVE_MEAN_WIND')
margins, at (`ABOVE_MEDIAN_RAIN')
