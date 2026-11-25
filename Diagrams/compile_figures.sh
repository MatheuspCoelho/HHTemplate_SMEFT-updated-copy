#!/usr/bin/env bash
#
# Compile every Feynman-diagram source in ./src into a cropped, standalone PDF
# stored in ./pdf. Run this once (and again whenever a diagram changes); the
# main document then simply \includegraphics the ready PDFs, so the Feynman
# diagrams are not re-rendered (and MetaPost is not needed) on every draft build.
#
# Usage: ./compile_figures.sh

set -euo pipefail

FIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${FIG_DIR}/.." && pwd)"
SRC_DIR="${FIG_DIR}/src"
PDF_DIR="${FIG_DIR}/pdf"
BUILD_DIR="${FIG_DIR}/build"

# The diagrams are compiled inside BUILD_DIR, so the project root has to be on
# the search path for the style files shipped with the paper (heppennames2.sty,
# hepparticles.sty) that provide the particle symbols used in the labels.
# ./texmf holds the standalone class (plus varwidth) that the diagram wrappers
# need and that the system TeX Live here does not provide.
export TEXINPUTS="${FIG_DIR}/texmf:${PROJECT_DIR}:${TEXINPUTS:-}"

# feynmp-auto needs MetaPost (mpost), which the system TeX Live does not ship.
# Fall back to the matching TeX Live 2020 mpost from CVMFS, exposed through a
# private bin directory holding just that one symlink. Putting the whole CVMFS
# bin directory on PATH would also override pdflatex and mix two TeX trees,
# which breaks the build with "Command \negmedspace already defined".
if ! command -v mpost >/dev/null 2>&1; then
  for mpdir in /cvmfs/sft.cern.ch/lcg/external/texlive/2020/bin/x86_64-linux; do
    if [[ -x "${mpdir}/mpost" ]]; then
      mkdir -p "${BUILD_DIR}/bin"
      ln -sf "${mpdir}/mpost" "${BUILD_DIR}/bin/mpost"
      export PATH="${BUILD_DIR}/bin:${PATH}"
      break
    fi
  done
fi

mkdir -p "${PDF_DIR}" "${BUILD_DIR}"

shopt -s nullglob
sources=("${SRC_DIR}"/*.tex)
if (( ${#sources[@]} == 0 )); then
  echo "No diagram sources found in ${SRC_DIR}" >&2
  exit 1
fi

for src in "${sources[@]}"; do
  name="$(basename "${src}" .tex)"
  wrapper="${BUILD_DIR}/${name}.tex"

  cat > "${wrapper}" <<TEXEOF
\documentclass[border={4pt 14pt},varwidth=17cm]{standalone}
% Deliberately NOT loading the paper's mathdesign/Charter fonts here: the
% MetaPost labels would then embed the same fonts the main document already
% uses, and the glyphs get substituted when the PDF is included, which shows up
% as garbled particle labels. Computer Modern keeps the included figures
% self-contained. amssymb supplies \Box for the SMEFT panel.
\usepackage{amssymb}
\usepackage{heppennames2}
\usepackage{feynmp-auto}
\usepackage{xcolor}
\usepackage{amsmath}
\usepackage{xspace}
% Category colours used by the diagrams, kept in sync with the definitions in
% the main document preamble.
\definecolor{cat1}{RGB}{87,144,252}
\definecolor{cat2}{RGB}{248,156,32}
\definecolor{cat3}{RGB}{228,37,54}
\definecolor{cat4}{RGB}{150,74,139}
\definecolor{cat6}{RGB}{122,33,221}
% Coupling macros (\klthree, ...) come from the paper's Definitions.
\input{${PROJECT_DIR}/Definitions}
\begin{document}
\input{${SRC_DIR}/${name}}
\end{document}
TEXEOF

  echo "==> compiling ${name}"
  # Build inside BUILD_DIR so feynmp-auto's mpost call and the generated .mp
  # files share the same working directory. Two passes: first writes/renders the
  # MetaPost figures, second includes them.
  ( cd "${BUILD_DIR}" \
      && pdflatex -shell-escape -interaction=nonstopmode -halt-on-error "${name}.tex" >/dev/null \
      && pdflatex -shell-escape -interaction=nonstopmode -halt-on-error "${name}.tex" >/dev/null )

  cp "${BUILD_DIR}/${name}.pdf" "${PDF_DIR}/${name}.pdf"
  echo "    -> pdf/${name}.pdf"
done

echo
echo "All diagrams compiled into ${PDF_DIR}"
