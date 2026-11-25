# `c_H` and `c_Hkin` dependence of single-Higgs rates

The figures of the section "Interplay with single-Higgs measurements" of the report.

| figure | made by |
|---|---|
| `c1_by_channel`, `kappa_lambda_vs_cH`, `production_ratio_vs_cH`, `br_ratio_vs_cH`, `sigma_times_br_vs_cH` | the SMEFT single-H model of multih-models, `SMEFT/singleH/scripts/plot_naive_cH.py` (recipe in `SMEFT/singleH/scripts/README.md` there) |
| `width_vs_cHkin`, `br_ratio_vs_cHkin`, `physical_region_cHkin`, `fccee_vs_lhc_vs_cH` | `single_higgs_ch.py`, this folder |

## Run

The multih-models script first: it writes the five figures above and
`single_h_naive_cH_scan.csv`, copied to `tables/`. Then

```bash
python3 single_higgs_ch.py                    # Lambda = 1 TeV
python3 single_higgs_ch.py --lambda-tev 2
python3 single_higgs_ch.py --repo-scan <file> # another scan of the multih-models model
```

It needs numpy and matplotlib. Outputs:

| output | content |
|---|---|
| `figures/width_vs_cHkin.pdf` | Γ_i/Γ_i,SM against `c_Hkin`, quadratic and linear truncation |
| `figures/br_ratio_vs_cHkin.pdf` | BR/BR_SM against `c_Hkin` |
| `figures/physical_region_cHkin.pdf` | the lower bound on `c_Hkin` for positive linearly truncated rates, against Λ |
| `figures/fccee_vs_lhc_vs_cH.pdf` | σ/σ_SM against `c_H`: `e+e- → ZH` at 240 GeV, ttH and ggH at 13 TeV |
| `tables/single_higgs_ch_<ecm>TeV.csv` | the standalone calculation below on a grid of `c_H` |

The ttH and ggH curves and the |κ_λ| > 20 bands of `fccee_vs_lhc_vs_cH` are read from
`tables/single_h_naive_cH_scan.csv`. The figures use the style of the multih-models ones (the
plain style of its `tools/plot_style.py`).

## What it computes

`c_H` enters through the naive translation onto the trilinear coupling (arXiv:2304.01968,
Tab. 2.1), `κ_λ = 1 − 2 (v²/Λ²)(v²/m_H²) c_H + 3 (v²/Λ²) c_Hkin`, and the NLO κ_λ dependence of
single-H rates (arXiv:1607.04251, 1709.08649), `Σ = Z_H Σ_LO (1 + κ_λ C₁ + δZ_H)`,
`Z_H = 1/(1 − (κ_λ² − 1) δZ_H)`, `δZ_H = −1.536 × 10⁻³`, normalised to its value at κ_λ = 1.
Branching ratios are the ratio of the partial widths to the total width; `Z_H` cancels in them.

`c_Hkin = c_HBox − c_HD/4` rescales every single-H coupling by `1 + x`, `x = c_Hkin v²/Λ²`, and
every partial width by `(1 + x)²`, `1 + 2x` at linear order; the linear form is zero at
`c_Hkin = −Λ²/(2v²)`, −8.25 at Λ = 1 TeV. The HH box-squared, interference and triangle-squared
terms scale as `(1 + x)⁴`, `(1 + x)⁶`, `(1 + x)⁸`; their linear forms are zero at
`c_Hkin = −Λ²/(p v²)`, p = 4, 6, 8 (`chkin_zero_crossing(Λ, p)`).

The `e+e-` production coefficients are in `C1_EE_PROD`: `C₁(ZH) = 1.5 %` at 240 GeV
(arXiv:1312.3322). `C₁(ZH)` at 365 GeV and `C₁(ννH)` are not filled (`None`); `c1_ee()` raises for
them and the figure omits their curves.

The production C₁ are tabulated at 7, 8, 13 and 14 TeV (arXiv:1607.04251, Tab. 2) and
interpolated linearly in between; above 14 TeV `c1_production()` raises unless `--allow-extrap`
(a flat extrapolation). The decay C₁ do not depend on √s.

## Inputs

| quantity | value | source |
|---|---|---|
| C₁ production, 13 TeV | ggF 0.66, VBF 0.64, WH 1.03, ZH 1.19, ttH 3.51 % | arXiv:1607.04251, Tab. 2 |
| C₁ decay | γγ 0.49, ZZ 0.83, WW 0.73, gg 0.66, ff̄ 0 % | arXiv:1607.04251, Tab. 1 |
| C₁(Zγ) | the γγ value (none in the reference); BR(Zγ) = 0.15 % | listed in `C1_DECAY_ASSUMED` |
| δZ_H | −1.536 × 10⁻³ | arXiv:1709.08649, Eq. (7) |
| range of the C₁ forms | \|κ_λ\| < 20 | arXiv:1709.08649 |
| SM branching ratios | m_H = 125.09 GeV | LHC Higgs XS WG, YR4 |
| v, m_H | 246.22 GeV, 125.25 GeV | |
| `c_H`, `c_Hkin` ranges | [−50, 50], [−8, 8] | the report, Tab. `SMEFTcouplings` |

The multih-models model differs in the production form (`Z_H (1 + (κ_λ − 1) C₁/K_EW)`), m_H
(125 GeV) and C₁(Zγ) (arXiv:1912.06429); its rates agree with the table of this script within
about 1 % over the `c_H` range.
