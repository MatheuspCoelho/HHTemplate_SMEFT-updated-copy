#!/usr/bin/env python3
"""
Self-coupling (c_H) dependence of single-Higgs rates via the naive HEFT->SMEFT
translation.

Reproduces the study of Sec. "Interplay with single-Higgs measurements" of the
HHTemplate_SMEFT report: the c_H dependence of single-Higgs production cross
sections, of the Higgs branching ratios, and of the product sigma x BR, over the
full c_H range used in the report.

Physics
-------
The trilinear coupling enters single-Higgs production and decay only through
electroweak loop corrections. Following Degrassi et al. (1607.04251) and
Maltoni et al. (1709.08649) these are captured by

    Sigma^BSM = Z_H^BSM * Sigma_LO * (1 + kappa_lambda * C1 + dZ_H)
    Z_H^BSM   = 1 / (1 - (kappa_lambda^2 - 1) * dZ_H)
    dZ_H      = -1.536e-3                                    (universal)

with C1 a process- and kinematics-dependent coefficient. The SMEFT Wilson
coefficients are mapped onto kappa_lambda with the naive translation of the
Warsaw-basis dictionary in Alasfar et al. (2304.01968, Table 2.1):

    kappa_lambda = 1 - 2 (v^2/Lambda^2)(v^2/mH^2) c_H + 3 (v^2/Lambda^2) c_Hkin

Limitations, stated in the report:
  * C1 assumes only the trilinear coupling is rescaled. A consistent SMEFT
    calculation at the same order has further one-loop pieces and operator
    mixing that this prescription does not capture.
  * The counting is inhomogeneous: an O(Lambda^-2) coefficient is inserted into
    a one-loop correction.
  * Eq. (9) of 1709.08649 is quoted as valid for |kappa_lambda| < 20.
    is_valid() flags points outside that.

Centre-of-mass energy
---------------------
C1 for production depends on sqrt(s) and is tabulated at 7, 8, 13 and 14 TeV
(1607.04251 Tab. 2); it decreases with energy.
  * 13 TeV     -> tabulated.
  * 13.6 TeV   -> interpolated linearly between 13 and 14 TeV; the 13 -> 14 TeV
                  change is at most 1.1 % relative (ttH).
  * 100 TeV    -> not computed. Requesting it raises unless allow_extrap (a
                  flat extrapolation).
C1 for decay does not depend on sqrt(s), nor do the branching ratios.

Usage
-----
    python3 single_higgs_ch.py                 # 13 TeV, writes figures + tables
    python3 single_higgs_ch.py --ecm 13.6
    python3 single_higgs_ch.py --lambda-tev 2

Outputs land in figures/ and tables/ next to this script.
"""

from __future__ import annotations

import argparse
import csv
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------
V_GEV = 246.22          # electroweak vev
MH_GEV = 125.25         # Higgs boson mass
DZ_H = -1.536e-3        # universal wave-function term, 1709.08649 Eq. (7)
KAPPA_VALIDITY = 20.0   # |kappa_lambda| range of validity, 1709.08649

# Production C1 [%] vs sqrt(s) [TeV], 1607.04251 Tab. 2.
# 1709.08649 Tab. 1 adds tHj and H->4l at 13 TeV only.
C1_PROD_VS_ECM = {
    #            7 TeV  8 TeV  13 TeV 14 TeV
    "ggF":      (0.66,  0.66,  0.66,  0.66),
    "VBF":      (0.65,  0.65,  0.64,  0.64),
    "WH":       (1.06,  1.05,  1.03,  1.03),
    "ZH":       (1.23,  1.22,  1.19,  1.18),
    "ttH":      (3.87,  3.78,  3.51,  3.47),
}
ECM_GRID = np.array([7.0, 8.0, 13.0, 14.0])
ECM_TABULATED_MAX = 14.0

# Decay C1 [%], 1607.04251 Tab. 1. Energy independent.
# ff-bar (bb, tautau, cc, mumu, ss) has C1 = 0 exactly at this order.
C1_DECAY = {
    "gamgam": 0.49,
    "ZZ": 0.83,
    "WW": 0.73,
    "gg": 0.66,
    "bb": 0.0,
    "tautau": 0.0,
    "cc": 0.0,
    "mumu": 0.0,
    "ss": 0.0,
    # Z-gamma is not given in 1607.04251 Tab. 1. Assigned the gamma-gamma value
    # as a stand-in; BR(Zgam) = 0.15% so this choice is numerically irrelevant
    # for the total width. Flagged rather than hidden.
    "Zgam": 0.49,
}
C1_DECAY_ASSUMED = {"Zgam"}

# --- e+e- production, for the FCC-ee comparison -----------------------------
# Only one value is sourced. McCullough (1312.3322) finds that the trilinear
# correction raises sigma(e+e- -> ZH) by about 1.5% for delta kappa_lambda = 1
# at sqrt(s) ~ 240-250 GeV, i.e. C1(ZH) ~ 1.5% in the convention used here.
#
# C1(ZH) at 365 GeV and C1(nunuH) are not filled (None): c1_ee() raises for
# them and fig_fccee() does not draw them. With a value from a reference, the
# curve is drawn.
C1_EE_PROD = {
    ("ZH", 240.0): 1.5,
    ("ZH", 365.0): None,     # TO FILL
    ("nunuH", 365.0): None,  # TO FILL
}


def c1_ee(channel: str, ecm_gev: float) -> float:
    """C1 [%] for an e+e- production channel at a given sqrt(s) in GeV."""
    try:
        val = C1_EE_PROD[(channel, float(ecm_gev))]
    except KeyError:
        raise KeyError(
            f"no C1 entry for {channel} at {ecm_gev} GeV; "
            f"known: {sorted(C1_EE_PROD)}"
        )
    if val is None:
        raise ValueError(
            f"C1({channel}) at {ecm_gev} GeV is not filled in -- see the note "
            "next to C1_EE_PROD. Refusing to guess."
        )
    return val


# SM branching ratios at mH = 125.09 GeV, LHC Higgs XS WG YR4.
BR_SM = {
    "bb": 0.5809,
    "WW": 0.2152,
    "gg": 0.0818,
    "tautau": 0.0627,
    "cc": 0.0289,
    "ZZ": 0.02641,
    "gamgam": 0.00227,
    "Zgam": 0.001541,
    "mumu": 0.0002171,
    "ss": 0.0002,
}

# Wilson-coefficient ranges used in the report (Tab. "SMEFTcouplings").
CH_RANGE = (-50.0, 50.0)
CHKIN_RANGE = (-8.0, 8.0)

# Petroff 6-colour scheme, as defined in the report preamble (used by the
# functions that are not run; the figures made use REPO_PALETTE).
CAT = ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#7a21dd", "#9c9ca1"]
INK = "#1a1a1a"
INK_MUTED = "#6b6b6b"
GRID = "#dcdcdc"


# --------------------------------------------------------------------------
# Physics
# --------------------------------------------------------------------------
def c1_production(channel: str, ecm_tev: float, allow_extrap: bool = False) -> float:
    """C1 [%] for a production channel at sqrt(s) = ecm_tev.

    Interpolates linearly on the published 7/8/13/14 TeV grid. Raises above
    14 TeV unless allow_extrap, because no calculation exists there.
    """
    if channel not in C1_PROD_VS_ECM:
        raise KeyError(f"no production C1 tabulated for {channel!r}")
    if ecm_tev > ECM_TABULATED_MAX and not allow_extrap:
        raise ValueError(
            f"C1 for production is published only up to {ECM_TABULATED_MAX} TeV "
            f"(1607.04251 Tab. 2); {ecm_tev} TeV would require a dedicated "
            f"calculation. allow_extrap=True gives a flat extrapolation; the "
            f"branching ratios do not depend on the energy."
        )
    return float(np.interp(ecm_tev, ECM_GRID, C1_PROD_VS_ECM[channel]))


def kappa_lambda(c_h, c_hkin=0.0, lambda_tev=1.0):
    """Naive translation, 2304.01968 Tab. 2.1 (Warsaw basis)."""
    lam_gev2 = (lambda_tev * 1000.0) ** 2
    v2_over_lam2 = V_GEV**2 / lam_gev2
    return (
        1.0
        - 2.0 * v2_over_lam2 * (V_GEV**2 / MH_GEV**2) * np.asarray(c_h, dtype=float)
        + 3.0 * v2_over_lam2 * np.asarray(c_hkin, dtype=float)
    )


def is_valid(kl):
    """Within the |kappa_lambda| < 20 validity of the C1 expansion."""
    return np.abs(kl) < KAPPA_VALIDITY


def z_h(kl):
    return 1.0 / (1.0 - (np.asarray(kl, dtype=float) ** 2 - 1.0) * DZ_H)


def rate_ratio(kl, c1_percent):
    """Sigma^BSM / Sigma^SM (or Gamma^BSM / Gamma^SM) for one channel."""
    c1 = c1_percent / 100.0
    kl = np.asarray(kl, dtype=float)
    num = z_h(kl) * (1.0 + kl * c1 + DZ_H)
    den = z_h(1.0) * (1.0 + c1 + DZ_H)
    return num / den


def br_ratio(kl, channel):
    """BR^BSM / BR^SM. Z_H cancels between numerator and total width, so the
    effect comes purely from the SPREAD in C1 across decay modes."""
    kl = np.asarray(kl, dtype=float)
    c1_i = C1_DECAY[channel] / 100.0
    partial = (1.0 + kl * c1_i + DZ_H) / (1.0 + c1_i + DZ_H)
    total = np.zeros_like(kl, dtype=float)
    for ch, br in BR_SM.items():
        c1_j = C1_DECAY[ch] / 100.0
        total = total + br * (1.0 + kl * c1_j + DZ_H) / (1.0 + c1_j + DZ_H)
    total = total / sum(BR_SM.values())
    return partial / total


def sigma_times_br(kl, prod, decay, ecm_tev, allow_extrap=False):
    return rate_ratio(kl, c1_production(prod, ecm_tev, allow_extrap)) * br_ratio(kl, decay)


# --------------------------------------------------------------------------
# c_Hkin = c_HBox - c_HD/4 : universal Higgs wave-function rescaling
# --------------------------------------------------------------------------
def width_universal_factor(c_hkin, lambda_tev=1.0, order="quadratic"):
    """Gamma_i / Gamma_i^SM from the c_Hkin field redefinition.

    Removing the c_HBox / c_HD contributions to the Higgs kinetic term rescales
    every single-Higgs coupling by (1 + x), x = c_Hkin v^2/Lambda^2, so every
    partial width scales by (1 + x)^2. Truncating that square at linear order
    leaves 1 + 2x, which turns NEGATIVE for x < -1/2, i.e.

        c_Hkin < -Lambda^2 / (2 v^2)          [see chkin_zero_crossing()]

    The quadratic form is a square and is not negative. The rescaling is the
    same for every partial width, so it cancels in the branching ratios.
    """
    x = np.asarray(c_hkin, dtype=float) * V_GEV**2 / (lambda_tev * 1000.0) ** 2
    if order == "linear":
        return 1.0 + 2.0 * x
    if order == "quadratic":
        return (1.0 + x) ** 2
    raise ValueError("order must be 'linear' or 'quadratic'")


def chkin_zero_crossing(lambda_tev=1.0, power=2):
    """c_Hkin at which a linearly truncated rate turns negative.

    An observable built from `power` factors of the universal coupling
    rescaling scales as (1+x)^power; truncated at linear order that is
    1 + power*x, which vanishes at x = -1/power. The relevant powers are

        power = 2  single-Higgs partial width      (c_i^2)
        power = 4  HH box-squared term             (c_t^4)
        power = 6  HH box-triangle interference    (c_t^3 c_hhh)
        power = 8  HH triangle-squared term        (c_t^2 c_hhh^2)

    since c_t = (1+x) and c_hhh = (1+x)^3: the HH bounds are 1/2 to 1/4 of the
    single-Higgs one.
    """
    return -((lambda_tev * 1000.0) ** 2) / (power * V_GEV**2)


# (label, power of (1+x)) for the observables that bound c_Hkin
CHKIN_TERMS = [
    (r"HH  $|A_{\Delta}|^2$", 8),
    (r"HH  interference", 6),
    (r"HH  $|A_{\Box}|^2$", 4),
    (r"single-$H$  $\Gamma_i$", 2),
]


# --------------------------------------------------------------------------
# Plot helpers
# --------------------------------------------------------------------------
def _style(ax, xlabel, ylabel, title=None):
    ax.set_xlabel(xlabel, fontsize=14, color=INK)
    ax.set_ylabel(ylabel, fontsize=14, color=INK)
    if title:
        ax.set_title(title, fontsize=13, color=INK, pad=8)
    ax.grid(True, color=GRID, lw=0.6, alpha=0.9)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK_MUTED)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK_MUTED, labelsize=12, length=3)


def _shade_invalid(ax, ch, kl):
    """Grey out where the C1 expansion is quoted as invalid."""
    bad = ~is_valid(kl)
    if not bad.any():
        return
    ymin, ymax = ax.get_ylim()
    ax.fill_between(ch, ymin, ymax, where=bad, color="#000000", alpha=0.055, lw=0)
    ax.set_ylim(ymin, ymax)
    # annotate once, on whichever invalid side exists
    for side in (ch[bad].min(), ch[bad].max()):
        if abs(side) == abs(ch[bad]).max():
            ax.text(
                side, ymax, r" $|\kappa_\lambda|>20$ ",
                fontsize=11, color=INK_MUTED, va="top",
                ha="left" if side < 0 else "right",
            )
            break


def _label_end(ax, x, y, text, color, dx=0.0, dy=0.0, ha="left"):
    ax.text(x + dx, y + dy, text, color=color, fontsize=12, va="center", ha=ha,
            fontweight="medium")


def _label_column(ax, x, items, min_gap_frac=0.062, ha="left"):  # noqa: D401
    """Direct-label several curves at the same x without letting the labels
    overlap: nudge them apart along y, keeping their original order.

    items: iterable of (y_data, text, color). ylim must already be final.
    """
    y0, y1 = ax.get_ylim()
    span = y1 - y0
    entries = sorted(
        ({"f": (y - y0) / span, "t": t, "c": c} for y, t, c in items),
        key=lambda e: e["f"],
    )
    # single forward pass pushing up, then clamp back down from the top
    for i in range(1, len(entries)):
        if entries[i]["f"] - entries[i - 1]["f"] < min_gap_frac:
            entries[i]["f"] = entries[i - 1]["f"] + min_gap_frac
    overflow = entries[-1]["f"] - 1.0
    if overflow > 0:
        for e in entries:
            e["f"] -= overflow
    for e in entries:
        ax.text(x, y0 + e["f"] * span, e["t"], color=e["c"], fontsize=12,
                va="center", ha=ha, fontweight="medium",
                bbox=dict(facecolor="white", edgecolor="none", pad=0.8, alpha=0.85))


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------
def _translation_note(ax, x=0.02, y=0.02, ha="left", va="bottom"):
    """Mark that the c_H dependence enters only through the naive translation."""
    ax.text(x, y, r"$c_H\!\to\!\kappa_\lambda$ naive transl.",
            transform=ax.transAxes, fontsize=8.5, color=INK_MUTED,
            ha=ha, va=va, zorder=6)


def fig_kappa(ch, outdir, lambdas=(1.0, 2.0, 3.0)):
    fig, ax = plt.subplots(figsize=(4.3, 3.3))
    curves = []
    for i, lam in enumerate(lambdas):
        kl = kappa_lambda(ch, 0.0, lam)
        ax.plot(ch, kl, color=CAT[i], lw=2.6, solid_capstyle="round")
        curves.append((kl[0], rf"$\Lambda={lam:g}$ TeV", CAT[i]))
    ax.axhline(1.0, color=INK_MUTED, lw=0.9, ls=(0, (4, 3)))
    _style(ax, r"$c_H$", r"$\kappa_\lambda$",
           r"$c_H \rightarrow \kappa_\lambda$ translation")
    ax.set_xlim(CH_RANGE[0], CH_RANGE[1])
    # grey out where the C1 expansion is not valid, as on the other panels
    y0, y1 = ax.get_ylim()
    for lo, hi in ((y0, -KAPPA_VALIDITY), (KAPPA_VALIDITY, y1)):
        if hi > lo:
            ax.axhspan(lo, hi, color="#000000", alpha=0.055, lw=0)
    ax.set_ylim(y0, y1)
    ax.text(CH_RANGE[1] * 0.97, KAPPA_VALIDITY, r"$|\kappa_\lambda|>20$",
            fontsize=11, color=INK_MUTED, va="bottom", ha="right")
    _label_column(ax, ch[0] + 0.03 * (ch[-1] - ch[0]), curves, min_gap_frac=0.085)
    fig.tight_layout()
    p = os.path.join(outdir, "kappa_lambda_vs_cH.pdf")
    fig.savefig(p); plt.close(fig)
    return p


def fig_production(ch, kl, ecm, outdir, allow_extrap=False):
    fig, ax = plt.subplots(figsize=(4.3, 3.3))
    # ggF and VBF have nearly identical C1 and so nearly identical curves;
    # they are given well-separated hues (purple / grey) rather than the two
    # adjacent purples, which would be the hardest pair to tell apart.
    order = [("ttH", CAT[0]), ("ZH", CAT[1]), ("WH", CAT[2]),
             ("ggF", CAT[3]), ("VBF", CAT[5])]
    labels = []
    for chan, col in order:
        c1 = c1_production(chan, ecm, allow_extrap)
        r = rate_ratio(kl, c1)
        ax.plot(ch, r, color=col, lw=2.6, solid_capstyle="round")
        labels.append((r[-1], f"  {chan}", col))
    ax.axhline(1.0, color=INK_MUTED, lw=0.9, ls=(0, (4, 3)))
    _style(ax, r"$c_H$", r"$\sigma/\sigma_{\rm SM}$",
           rf"Production, $\sqrt{{s}}={ecm:g}$ TeV")
    ax.set_xlim(CH_RANGE[0], CH_RANGE[1] * 1.24)
    _shade_invalid(ax, ch, kl)
    _label_column(ax, ch[-1], labels)
    _translation_note(ax)
    fig.tight_layout()
    p = os.path.join(outdir, "production_ratio_vs_cH.pdf")
    fig.savefig(p); plt.close(fig)
    return p


def read_repo_scan(path):
    """The cH scan of the multih-models SMEFT single-H model
    (single_h_naive_cH_scan.csv of SMEFT/singleH/scripts/plot_naive_cH.py):
    {column: array}, columns cH, kl and xs_<mode>, Lambda = 1 TeV, 13 TeV."""
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def fig_fccee(outdir, repo_scan, ecm_pp=13.0):
    """sigma/sigma_SM vs c_H for e+e- -> ZH at FCC-ee, against pp channels.

    The pp curves (ttH, ggH) and kappa_lambda are those of the multih-models
    model (repo_scan); the e+e- ones use the C1 of C1_EE_PROD in rate_ratio.
    The e+e- channels without a C1 entry are not drawn.
    """
    _repo_style()
    scan = read_repo_scan(repo_scan)
    ch, kl = scan["cH"], scan["kl"]
    fig, ax = plt.subplots(figsize=(5, 4))
    i = 0
    for (chan, ecm), tex in ((("ZH", 240.0), r"$ZH$"), (("ZH", 365.0), r"$ZH$"),
                             (("nunuH", 365.0), r"$\nu\bar{\nu}H$")):
        try:
            c1 = c1_ee(chan, ecm)
        except (KeyError, ValueError):
            continue
        ax.plot(ch, rate_ratio(kl, c1), color=REPO_PALETTE[1 + i], lw=2.2,
                label=r"$e^+e^-\to$" + tex + ", {:g} GeV".format(ecm))
        i += 1
    for mode, color, tex in (("ttH", REPO_PALETTE[0], r"$t\bar{t}H$"),
                             ("ggH", REPO_PALETTE[3], "ggH")):
        ax.plot(ch, scan["xs_" + mode], color=color, lw=2.2, ls="--",
                label=tex + ", {:g} TeV".format(ecm_pp))
    ax.axhline(1.0, color="0.4", lw=0.8, ls=":")
    ax.set_xlim(CH_RANGE[0], CH_RANGE[1])
    y0, y1 = ax.get_ylim()
    ax.fill_between(ch, y0, y1, where=~is_valid(kl), color="k", alpha=0.06, lw=0)
    ax.set_ylim(y0, y1)
    ax.set_xlabel(r"$c_{H}$")
    ax.set_ylabel(r"$\sigma/\sigma_{\mathrm{SM}}$")
    ax.legend(frameon=False, fontsize="small", loc="lower left")
    fig.tight_layout()
    p = os.path.join(outdir, "fccee_vs_lhc_vs_cH.pdf")
    fig.savefig(p); plt.close(fig)
    return p


def fig_br(ch, kl, outdir):
    fig, ax = plt.subplots(figsize=(4.3, 3.3))
    order = ["ZZ", "WW", "gg", "gamgam", "bb"]
    tex = {"ZZ": r"$ZZ$", "WW": r"$WW$", "gg": r"$gg$",
           "gamgam": r"$\gamma\gamma$", "bb": r"$b\bar b,\ \tau\tau$"}
    labels = []
    for i, d in enumerate(order):
        r = br_ratio(kl, d)
        ax.plot(ch, r, color=CAT[i], lw=2.6, solid_capstyle="round")
        labels.append((r[-1], "  " + tex[d], CAT[i]))
    ax.axhline(1.0, color=INK_MUTED, lw=0.9, ls=(0, (4, 3)))
    _style(ax, r"$c_H$", r"${\rm BR}/{\rm BR}_{\rm SM}$",
           r"Branching ratios ($\sqrt{s}$ independent)")
    ax.set_xlim(CH_RANGE[0], CH_RANGE[1] * 1.30)
    _shade_invalid(ax, ch, kl)
    _label_column(ax, ch[-1], labels, min_gap_frac=0.082)
    _translation_note(ax)
    fig.tight_layout()
    p = os.path.join(outdir, "br_ratio_vs_cH.pdf")
    fig.savefig(p); plt.close(fig)
    return p


def fig_sigma_br(ch, kl, ecm, outdir, allow_extrap=False):
    fig, ax = plt.subplots(figsize=(4.3, 3.3))
    combos = [
        ("ttH", "gamgam", r"$t\bar tH(\gamma\gamma)$"),
        ("ttH", "bb", r"$t\bar tH(b\bar b)$"),
        ("ttH", "ZZ", r"$t\bar tH(ZZ)$"),
        ("ggF", "gamgam", r"$ggF(\gamma\gamma)$"),
    ]
    labels = []
    for i, (prod, dec, lab) in enumerate(combos):
        r = sigma_times_br(kl, prod, dec, ecm, allow_extrap)
        ax.plot(ch, r, color=CAT[i], lw=2.6, solid_capstyle="round")
        labels.append((r[-1], "  " + lab, CAT[i]))
    ax.axhline(1.0, color=INK_MUTED, lw=0.9, ls=(0, (4, 3)))
    _style(ax, r"$c_H$", r"$(\sigma\times{\rm BR})/{\rm SM}$",
           rf"Observable rate, $\sqrt{{s}}={ecm:g}$ TeV")
    ax.set_xlim(CH_RANGE[0], CH_RANGE[1] * 1.52)
    _shade_invalid(ax, ch, kl)
    _label_column(ax, ch[-1], labels)
    _translation_note(ax)
    fig.tight_layout()
    p = os.path.join(outdir, "sigma_times_br_vs_cH.pdf")
    fig.savefig(p); plt.close(fig)
    return p


# The style of the figures made with the multih-models SMEFT single-H model
# (SMEFT/singleH/scripts/plot_naive_cH.py there): matplotlib default, the plain
# style of its tools/plot_style.py, a legend and no title.
REPO_PALETTE = ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"]


def _repo_style():
    plt.style.use("default")
    plt.rcParams.update({"font.size": 13, "axes.grid": True, "grid.alpha": 0.3,
                         "axes.prop_cycle": plt.cycler(color=REPO_PALETTE)})


def fig_width_chkin(outdir, lambdas=(1.0,), npoints=1201):
    """Partial widths vs c_Hkin: where the linear truncation goes negative."""
    _repo_style()
    lam = lambdas[0]
    x0 = chkin_zero_crossing(lam)
    xlo, xhi = min(CHKIN_RANGE[0], x0) - 1.5, CHKIN_RANGE[1]
    chk = np.linspace(xlo, xhi, npoints)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(chk, width_universal_factor(chk, lam, "quadratic"), color=REPO_PALETTE[0], lw=2.2,
            label="quadratic")
    ax.plot(chk, width_universal_factor(chk, lam, "linear"), color=REPO_PALETTE[2], lw=2.2,
            ls="--", label="linear")
    ax.axhline(1.0, color="0.4", lw=0.8, ls="--")
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xlim(xlo, xhi)
    y0, y1 = ax.get_ylim()
    y0 = min(y0, -0.35)
    ax.axhspan(y0, 0.0, color="k", alpha=0.06, lw=0)
    ax.set_ylim(y0, y1)
    ax.axvline(CHKIN_RANGE[0], color="0.4", lw=0.8, ls=":")
    ax.plot([x0], [0.0], "o", color=REPO_PALETTE[2], ms=6)
    ax.annotate(rf"${x0:.2f}$", xy=(x0, 0.0), xytext=(x0 + 0.6, 0.35), fontsize="small",
                color=REPO_PALETTE[2])
    ax.set_xlabel(r"$c_{H,\mathrm{kin}}$")
    ax.set_ylabel(r"$\Gamma_i/\Gamma_i^{\mathrm{SM}}$")
    ax.legend(frameon=False, fontsize="small", loc="upper left")
    fig.tight_layout()
    p = os.path.join(outdir, "width_vs_cHkin.pdf")
    fig.savefig(p); plt.close(fig)
    return p


def fig_br_chkin(outdir, lambda_tev=1.0, npoints=801):
    """BR vs c_Hkin. The universal rescaling cancels between the partial and the
    total widths; the BR change through the kappa_lambda term of the translation."""
    _repo_style()
    chk = np.linspace(CHKIN_RANGE[0], CHKIN_RANGE[1], npoints)
    kl = kappa_lambda(0.0, chk, lambda_tev)
    tex = {"ZZ": r"$ZZ$", "WW": r"$WW$", "gg": r"$gg$", "gamgam": r"$\gamma\gamma$",
           "bb": r"$b\bar{b}$, $\tau\tau$"}
    fig, ax = plt.subplots(figsize=(5, 4))
    for i, d in enumerate(["ZZ", "WW", "gg", "gamgam", "bb"]):
        ax.plot(chk, br_ratio(kl, d), color=REPO_PALETTE[i], lw=2.2, label=tex[d])
    ax.axhline(1.0, color="0.4", lw=0.8, ls="--")
    ax.set_xlim(*CHKIN_RANGE)
    ax.set_xlabel(r"$c_{H,\mathrm{kin}}$")
    ax.set_ylabel(r"$\mathrm{BR}/\mathrm{BR}_{\mathrm{SM}}$")
    ax.legend(frameon=False, fontsize="small")
    fig.tight_layout()
    p = os.path.join(outdir, "br_ratio_vs_cHkin.pdf")
    fig.savefig(p); plt.close(fig)
    return p


def fig_physical_region(outdir, lam_range=(0.5, 2.0), npoints=601):
    """The lower bound on c_Hkin for positive linearly truncated rates, against Lambda.

    With the linear truncation Gamma_i/Gamma_i^SM = 1 + 2x, so the boundary is
    c_Hkin > -Lambda^2 / (2 v^2); the strictest HH term (triangle squared) gives
    c_Hkin > -Lambda^2 / (8 v^2).
    """
    _repo_style()
    lam = np.linspace(*lam_range, npoints)
    floor_br = -(lam * 1000.0) ** 2 / (2.0 * V_GEV**2)
    floor_hh = -(lam * 1000.0) ** 2 / (8.0 * V_GEV**2)
    ylo, yhi = -20.0, 12.0
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.fill_between(lam, ylo, floor_br, color="k", alpha=0.06, lw=0)
    ax.plot(lam, floor_br, color=REPO_PALETTE[2], lw=2.2,
            label=r"single $H$: $\Gamma_i$, BR $> 0$")
    ax.plot(lam, floor_hh, color=REPO_PALETTE[5], lw=2.2, ls="--",
            label=r"HH: $|A_\Delta|^2 > 0$")
    for y in CHKIN_RANGE:
        ax.axhline(y, color="0.4", lw=0.8, ls=":")
    for power, color in ((2, REPO_PALETTE[2]), (8, REPO_PALETTE[5])):
        y = chkin_zero_crossing(1.0, power)
        ax.plot([1.0], [y], "o", color=color, ms=6)
        ax.annotate(rf"${y:.2f}$", xy=(1.0, y), xytext=(1.04, y + 0.8), fontsize="small",
                    color=color)
    ax.set_xlim(*lam_range)
    ax.set_ylim(ylo, yhi)
    ax.set_xlabel(r"$\Lambda$ [TeV]")
    ax.set_ylabel(r"$c_{H,\mathrm{kin}}$")
    ax.legend(frameon=False, fontsize="small", loc="lower left")
    fig.tight_layout()
    p = os.path.join(outdir, "physical_region_cHkin.pdf")
    fig.savefig(p); plt.close(fig)
    return p


def fig_c1_by_channel(outdir, ecm=13.0):
    """C1 per channel -- replaces the tabular listing."""
    prod = [(c, c1_production(c, ecm)) for c in ["ttH", "ZH", "WH", "ggF", "VBF"]]
    prod.append(("tHj", 0.91))
    dec = [(r"$H\to ZZ$", 0.83), (r"$H\to WW$", 0.73), (r"$H\to gg$", 0.66),
           (r"$H\to\gamma\gamma$", 0.49), (r"$H\to f\bar f$", 0.0)]
    prod.sort(key=lambda t: t[1])
    dec.sort(key=lambda t: t[1])

    fig, ax = plt.subplots(figsize=(4.6, 3.5))
    ypos, labels = [], []
    y = 0
    for group, col in ((dec, CAT[3]), (prod, CAT[0])):
        for name, val in group:
            ax.barh(y, val, height=0.6, color=col, zorder=3)
            ax.text(val + 0.10, y, f"{val:.2f}", va="center", ha="left",
                    fontsize=11, color=INK)
            ypos.append(y); labels.append(name)
            y += 1
        y += 0.7

    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=11.5)
    _style(ax, r"$C_1$  [%]", "", rf"$C_1$ by channel, $\sqrt{{s}}={ecm:g}$ TeV")
    ax.set_xlim(0, 4.4)
    ax.set_ylim(-0.8, y - 0.3)
    ax.grid(axis="y", visible=False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=CAT[0]),
               plt.Rectangle((0, 0), 1, 1, color=CAT[3])]
    ax.legend(handles, ["production", "decay"], loc="lower right", fontsize=11,
              frameon=False)
    fig.tight_layout(pad=0.6)
    p = os.path.join(outdir, "c1_by_channel.pdf")
    fig.savefig(p); plt.close(fig)
    return p


def write_table(ecm, outdir, lambda_tev, allow_extrap=False):
    rows = []
    for c in [-50, -40, -30, -20, -10, -5, -3, -2, -1, 0, 1, 2, 3, 5, 10, 20, 30, 40, 50]:
        kl = kappa_lambda(c, 0.0, lambda_tev)
        row = {
            "c_H": c,
            "kappa_lambda": round(float(kl), 4),
            "valid_kappa_lt_20": bool(is_valid(kl)),
        }
        for chan in ["ggF", "VBF", "WH", "ZH", "ttH"]:
            row[f"sigma_{chan}"] = round(float(rate_ratio(kl, c1_production(chan, ecm, allow_extrap))), 5)
        for d in ["gamgam", "ZZ", "WW", "bb", "gg"]:
            row[f"BR_{d}"] = round(float(br_ratio(kl, d)), 5)
        row["ttH_x_gamgam"] = round(float(sigma_times_br(kl, "ttH", "gamgam", ecm, allow_extrap)), 5)
        row["ttH_x_bb"] = round(float(sigma_times_br(kl, "ttH", "bb", ecm, allow_extrap)), 5)
        rows.append(row)
    p = os.path.join(outdir, f"single_higgs_ch_{ecm:g}TeV.csv")
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return p, rows


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ecm", type=float, default=13.0,
                    help="centre-of-mass energy in TeV (production C1 known to 14)")
    ap.add_argument("--lambda-tev", type=float, default=1.0, help="EFT scale in TeV")
    ap.add_argument("--allow-extrap", action="store_true",
                    help="production C1 above 14 TeV as a flat extrapolation")
    ap.add_argument("--npoints", type=int, default=801)
    ap.add_argument("--repo-scan", default=os.path.join(HERE, "tables", "single_h_naive_cH_scan.csv"),
                    help="the cH scan of the multih-models SMEFT single-H model, for the pp "
                         "curves of the FCC-ee figure (README.md)")
    args = ap.parse_args()

    figdir = os.path.join(HERE, "figures")
    tabdir = os.path.join(HERE, "tables")
    os.makedirs(figdir, exist_ok=True)
    os.makedirs(tabdir, exist_ok=True)

    ch = np.linspace(CH_RANGE[0], CH_RANGE[1], args.npoints)
    kl = kappa_lambda(ch, 0.0, args.lambda_tev)

    print(f"sqrt(s) = {args.ecm} TeV, Lambda = {args.lambda_tev} TeV, "
          f"c_H in [{CH_RANGE[0]:g}, {CH_RANGE[1]:g}]")
    print(f"  kappa_lambda spans [{kl.min():.2f}, {kl.max():.2f}]")
    frac = 100.0 * is_valid(kl).mean()
    print(f"  {frac:.0f}% of the c_H range satisfies |kappa_lambda| < 20")
    if C1_DECAY_ASSUMED:
        print(f"  note: C1 assumed (not published) for {sorted(C1_DECAY_ASSUMED)}")
    print(f"  note: C1(production) is sqrt(s) dependent, published 7-14 TeV; "
          f"C1(decay) and hence all BR results are sqrt(s) independent.")
    x0 = chkin_zero_crossing(args.lambda_tev)
    edge = width_universal_factor(CHKIN_RANGE[0], args.lambda_tev, "linear")
    edge_q = width_universal_factor(CHKIN_RANGE[0], args.lambda_tev, "quadratic")
    print(f"  c_Hkin: linearly truncated partial widths vanish at "
          f"c_Hkin = {x0:.2f} (Lambda = {args.lambda_tev:g} TeV)")
    print(f"          at the report's lower edge c_Hkin = {CHKIN_RANGE[0]:g}: "
          f"Gamma/Gamma_SM = {edge:.3f} (linear) vs {edge_q:.3f} (quadratic)")
    if x0 > CHKIN_RANGE[0]:
        print(f"          WARNING: zero crossing is INSIDE the scanned range")
    print(f"  physical range required of c_Hkin (BR/widths positive, linear "
          f"truncation): c_Hkin > {chkin_zero_crossing(args.lambda_tev, 2):.2f} "
          f"at Lambda = {args.lambda_tev:g} TeV")

    # c1_by_channel, kappa_lambda_vs_cH, production_ratio_vs_cH, br_ratio_vs_cH and
    # sigma_times_br_vs_cH of the report are made with the multih-models SMEFT single-H
    # model (README.md); their functions here are kept and not run
    outs = [
        fig_fccee(figdir, args.repo_scan),
        fig_width_chkin(figdir, lambdas=(args.lambda_tev,)),
        fig_br_chkin(figdir, args.lambda_tev),
        fig_physical_region(figdir),
    ]
    tab, rows = write_table(args.ecm, tabdir, args.lambda_tev, args.allow_extrap)
    outs.append(tab)
    for o in outs:
        print("  wrote", os.path.relpath(o, HERE))

    print("\n  c_H   kappa_l   sig(ttH)   BR(gamgam)   ttH x gamgam")
    for r in rows:
        if r["c_H"] in (-50, -10, -1, 0, 1, 10, 50):
            print(f"  {r['c_H']:>4}  {r['kappa_lambda']:>7.2f}   "
                  f"{r['sigma_ttH']:>8.4f}   {r['BR_gamgam']:>10.4f}   "
                  f"{r['ttH_x_gamgam']:>12.4f}")


if __name__ == "__main__":
    main()
