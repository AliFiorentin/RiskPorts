# CLAUDE.md

Guidance for Claude Code and other AI assistants working in this repository.

## What this repository is

Replication material for a *Marine Policy* submission estimating the
opportunity cost of weather-related stoppages at the São Francisco do Sul port
complex (Brazil). It is a **research archive**, not an evolving software
project: the code exists to reproduce published numbers, not to be refactored.

## The single most important rule

**Do not change the estimation.** Formulas, matching specifications,
estimands and variable sets in `scripts/01`–`03` are the ones that produced the
published tables. Improving them silently makes the archive worthless. If
something looks statistically wrong, say so — do not fix it in place.

Presentation is a different matter: figure styling, labels, resolution and
documentation can be improved freely. That is what `scripts/04` is for.

## Data

`data/baserisk1.dta` is the **only** analysis base. 4,332 vessel calls,
2018–2021. It reproduces Tables 1, 2 and A.1 of the paper exactly — verified
coefficient by coefficient.

The upstream source cannot be re-downloaded: ANTAQ restructured its statistical
system and dropped the revenue variable. There is no way to rebuild this file.
Treat it as irreplaceable.

Read `data/README.md` before touching any variable. It documents the cargo
dummy mapping (which does not follow the obvious 1–4 order), the missingness
pattern that drives sample size, and known data quality issues.

## Things that will trip you up

**The cargo dummies are not in listed order.** `d_carga1` is containerised
cargo, `d_carga2` general cargo, `d_carga3` liquid bulk, `d_carga4` dry bulk.
Verified by cross-tabulation, and it contradicts how the categories are ordered
in the paper's text.

**Weather threshold dummies use hard-coded cut-offs.** `d1*` uses the sample
median and `d2*` the sample mean, both written as numeric literals. Change the
sample and the cut-offs silently become wrong.

**Matching with replacement is not deterministic.** Ties are broken at random
and the original scripts set no seed. Point estimates reproduce exactly;
Rosenbaum bounds do not. Set a seed before comparing runs.

**`love.plot` with two statistics returns a `gtable`, not a `ggplot`.** Adding
`+ theme(...)` to it produces a blank figure with no error. Apply themes with
`theme_set()` beforehand.

**Times are in hours.** If a time variable shows values in the thousands, the
file is a corrupted export — stop and check.

## Known defects in the original scripts

Left in place deliberately, with fixes limited to what prevented execution:

- `01` mixes `MatchIt` and `Matching`; `base_pareada` is assigned twice, the
  second time from a `Matching::Match` object.
- `02` labels an `optimal, ratio = 2` call as "radius matching".
- `03` has no `use` statement — it operates on whatever is already in memory.
- The paper's text claims a caliper of 0.001 and clustered standard errors;
  the code applies neither. Standard errors in Table 2 are ordinary OLS.

## Conventions

Comments and messages are in English. Scripts are ASCII-only, to avoid encoding
problems between Windows and Linux. Paths are resolved relative to the
repository root — do not add `setwd()`.
