# Codebook — `baserisk1.dta`

Stata 118 format. **4,332 rows** (one per vessel call) × **79 columns**.
Coverage: 2018–2021, São Francisco do Sul port complex, Santa Catarina, Brazil.

Only 18 columns are source data; the remaining 61 are derived inside
`scripts/01_matching_psm.R` (squared terms and threshold dummies).

## Identification and facility

| Column | Type | Description |
|---|---|---|
| `ID` | num | Vessel call identifier (1,111–5,442) |
| `ano`, `mês` | num, text | Year and month |
| `chegada`, `desatrac` | text | Arrival and unberthing dates, `dd/mm/yyyy` |
| `instalacao` | text | Port facility |
| `t_instalacao` | text | `Porto Público` / `Porto Privado` |
| `terminal` | text | Terminal / berth |
| `imo` | num | Vessel IMO number (23 missing) |

Facilities: Porto Itapoá Terminais Portuários (2,129 calls), São Francisco do
Sul — Cais Público (1,115) and TESC (696), Terminal Aquaviário de São Francisco
do Sul (392).

## Treatment and outcomes

| Column | Description |
|---|---|
| `parada` | **Treatment.** 1 if the call experienced a weather-related stoppage (1,280 calls, 29.5%) |
| `t_atrac` | Waiting time to berth (hours) |
| `t_inicioop` | Waiting time to start operations (hours) |
| `t_operac` | Operation time (hours) |
| `t_desatrac` | Unberthing time (hours) |
| `t_atracado` | Time at berth (hours) |
| `t_estadia` | Total port stay (hours) |

## Cargo, port and revenue

| Column | Description |
|---|---|
| `movimentacao` | Cargo volume (tonnes), 0.03–281,258 |
| `perfil_carga` | Cargo profile, text |
| `p_carga` | Cargo profile, coded 1–4 |
| `d_carga1`–`d_carga4` | Cargo dummies — see mapping below |
| `terminais` | Terminal, coded 1–4 |
| `d_terminal1`–`d_terminal4` | Terminal dummies |
| `p_publico` | 1 if public port |
| `receita` | Revenue (BRL). **Missing in 2,521 of 4,332 rows** |
| `p_horas` | Auxiliary stoppage indicator |

### Cargo dummy mapping

Verified by cross-tabulation against `perfil_carga` (exact 1:1 correspondence).
Note that the mapping does **not** follow the order in which the categories are
listed in the paper:

| Dummy | Cargo profile | n |
|---|---|---|
| `d_carga1` | Containerised cargo | 2,130 |
| `d_carga2` | General cargo | 852 |
| `d_carga3` | Liquid and gas bulk | 399 |
| `d_carga4` | Dry bulk | 951 |

## Weather and waves

Matched to each vessel call. `mm_*` is rainfall (mm), `umid_*` relative
humidity (%), `vmax_*` maximum wind speed, for the Joinville, Barra Velha
(`barra`) and Itapoá stations. `mm_itapoa_e` is an additional Itapoá rainfall
series. `hs_medEst` is significant wave height at the harbour entrance.

Missingness varies: `mm_*` and `vmax_barra` are complete (4,332);
`umid_joiville` and `vmax_joiville` have 4,312; `hs_medEst` 4,320;
`umid_barra` 4,126; `umid_itapoa` 3,155; `vmax_itapoa` 2,671.

`vmax_itapoa` is the binding constraint on sample size: the regressions in
Table 2 of the paper run on 1,001 observations, down from 1,834 in the matched
sample, mostly because of it.

## Derived columns

`movimentacao2` and the `*2` suffixes are squared terms. `d_ano1`–`d_ano5` are
year dummies. `d1*` and `d2*` are weather threshold dummies: `d1` marks values
at or above the sample **median**, `d2` at or above the sample **mean**. Both
thresholds are hard-coded as literals in `scripts/01_matching_psm.R` — they are
not recomputed if the sample changes.

## Known data quality issues

These predate this work and are documented for transparency.

- `d_ano1` is constant zero across all 4,332 rows, yet enters the regression
  specifications. R drops it automatically, so estimates are unaffected.
- Two rows carry implausible arrival dates (03/01/2008 and 26/05/2015).
- Rainfall maxima are implausibly high for daily accumulation
  (`mm_joiville` reaches 12,491 mm) and disagree with the maxima reported in
  Table 1 of the paper, although means and standard deviations match exactly.
- `vmax_joiville` reaches 142.59 m/s, which is not physically plausible.

## Provenance

An earlier export of this dataset circulated with a corrupted decimal
separator: roughly 29% of rows had the decimal point dropped from the time
variables, so values such as `23.1667` hours appeared as `231667`. That file is
**not** included here. `baserisk1.dta` is the intact version and reproduces
Tables 1, 2 and A.1 of the paper exactly.
