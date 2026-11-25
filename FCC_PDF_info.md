PDF for FCC hh

Sources:
- PDF4LHC21 combination (arXiv:2203.05506) (https://arxiv.org/abs/2203.05506)
- Physics at a 100 TeV pp collider: Higgs and EW (arXiv:1606.09408) (https://arxiv.org/pdf/1606.09408)
- Higgs self-coupling measurements at the FCC-hh (arXiv:2312.03513) (https://arxiv.org/pdf/2312.03513)
- Measuring the Higgs self-coupling via HH at 100 TeV (EPJC) (https://link.springer.com/article/10.1140/epjc/s10052-020-08595-3)
- NNPDF3.1sx — small-x resummation (https://nnpdf.mi.infn.it/nnpdf3-1sx/)
- The small-x gluon for a 100 TeV collider (arXiv:1705.04217) (https://arxiv.org/pdf/1705.04217)

For FCC-hh the order-matched central choice (the POWHEG ME is NLO) is NNPDF4.0 NLO. 

1. Small-x gluon. At 13.6 TeV, gg→HH sits around x ∼ 0.01–0.1 where the gluon is well measured. At 100 TeV the same final state pulls the dominant sampling down to x ∼ few×10⁻³ and into the 10⁻⁵ tail — exactly where NNPDF4.0's gluon is weakly constrained by data (largely extrapolation below x∼10⁻³). So the central value and especially the PDF uncertainty degrade relative to 13.6 TeV. The literature is explicit that PDFs become the dominant theory uncertainty for EW-scale production at 100 TeV, precisely from this small-x region.
2. What FCC-hh actually uses. FCC-hh Higgs / di-Higgs self-coupling studies conventionally run NNPDF3.1/4.0 NLO for the central prediction and take PDF4LHC21 (CT18⊕MSHT20⊕NNPDF3.1) for the uncertainty band. So NNPDF4.0 is in the right family — you're not using something exotic — but the community estimates the uncertainty with the PDF4LHC21 combination, not a single set.
3. Small-x resummation. At FCC energies, NLLx small-x resummation (NNPDF3.1sx) shifts the small-x gluon noticeably vs the fixed-order fit. Not needed for a first look, but relevant before quoting final 100 TeV numbers.

Recommendation

- For shapes / a first cross-section estimate: keeping 90400 (NNPDF4.0 NNLO) is defensible — I would not block on it. Shapes are far more robust than the normalization.
- Before trusting absolute cross-sections at 100 TeV: (a) propagate the full replica PDF uncertainty (expect it to be visibly larger than at 13.6 TeV), (b) cross-check the central value/uncertainty against PDF4LHC21, and (c) consider a NNPDF3.1sx comparison for the small-x tail.
- One convention nit: the ME here is NLO (POWHEG), so if you want to match orders you'd use an NLO PDF (NNPDF4.0 NLO) rather than NNLO — though CMS Run-3 convention pairs NNLO PDFs with these, so 90400 is consistent with that convention if you're staying CMS-aligned.


PDF per ECM (now automatic)

┌────────────────┬────────────┬───────────┬─────────────┬────────────────────────────────────────────────────────┐
│     --ecm      │ ebeam/beam │ comEnergy │ PDF (LHAID) │                        PDF name                        │
├────────────────┼────────────┼───────────┼─────────────┼────────────────────────────────────────────────────────┤
│ 13             │ 6500       │ 13000     │ 90400       │ PDF4LHC15_nlo_30_pdfas (template default, unchanged)   │
├────────────────┼────────────┼───────────┼─────────────┼────────────────────────────────────────────────────────┤
│ 13.6 (default) │ 6800       │ 13600     │ 90400       │ PDF4LHC15_nlo_30_pdfas (current production, unchanged) │
├────────────────┼────────────┼───────────┼─────────────┼────────────────────────────────────────────────────────┤
│ 100 (FCC-hh)   │ 50000      │ 100000    │ 93300       │ PDF4LHC21_40_pdfas (auto-selected)                     │
└────────────────┴────────────┴───────────┴─────────────┴────────────────────────────────────────────────────────┘

Why 93300 for FCC-hh: it's the modern PDF4LHC combination (direct successor to the PDF4LHC15 the workflow already uses), NLO — matching the NLO POWHEG ME, the _pdfas variant like the current set, and its grid (XMin 1.0×10⁻⁶, QMax 10⁶ GeV) safely covers FCC-hh small-x/high-Q — whereas PDF4LHC15's XMin 6.1×10⁻⁶ would hit the grid edge at 100 TeV. The PDF is baked into the card's lhans1/lhans2 at gridpack-build time (written as a bare integer), so nanogen inherits it automatically — no PDF flag needed downstream. Override any energy with submit_smeft.sh --pdf <LHAID>.

Still recommended before trusting absolute 100 TeV cross-sections: quote the full PDF-replica uncertainty (large at small-x) and consider an NNPDF3.1sx small-x-resummation cross-check.