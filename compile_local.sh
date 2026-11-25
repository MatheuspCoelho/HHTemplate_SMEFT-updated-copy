#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

MAIN_TEX="${1:-SciPostPhysCommRep.tex}"
OUT_DIR="${2:-build-local}"

if [[ "${MAIN_TEX}" == "-h" || "${MAIN_TEX}" == "--help" ]]; then
  cat <<'EOF'
Usage: ./compile_local.sh [main-tex-file] [output-dir]

Defaults:
  main-tex-file: SciPostPhysCommRep.tex
  output-dir:    build-local

Examples:
  ./compile_local.sh
  ./compile_local.sh SciPostPhysCommRep.tex build-local
EOF
  exit 0
fi

if [[ ! -f "${MAIN_TEX}" ]]; then
  echo "Error: TeX source file not found: ${MAIN_TEX}" >&2
  exit 1
fi

for cmd in latexmk pdflatex bibtex kpsewhich; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Error: required command not found: ${cmd}" >&2
    exit 1
  fi
done

required_packages=(
  tocloft.sty
  tocbibind.sty
  titlesec.sty
  doi.sty
  lineno.sty
  mathdesign.sty
  feynmp-auto.sty
  geometry.sty
  hyperref.sty
  caption.sty
  cite.sty
)

missing_packages=()
for package in "${required_packages[@]}"; do
  if ! kpsewhich "${package}" >/dev/null 2>&1; then
    missing_packages+=("${package}")
  fi
done

if (( ${#missing_packages[@]} > 0 )); then
  echo "Error: missing required LaTeX packages:" >&2
  printf '  - %s\n' "${missing_packages[@]}" >&2
  echo "Install the missing TeX Live packages locally, then rerun the script." >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"

latexmk -C -output-directory="${OUT_DIR}" "${MAIN_TEX}" >/dev/null 2>&1 || true
latexmk \
  -pdf \
  -shell-escape \
  -interaction=nonstopmode \
  -file-line-error \
  -output-directory="${OUT_DIR}" \
  "${MAIN_TEX}"

pdf_name="$(basename "${MAIN_TEX%.tex}.pdf")"
echo
echo "Build completed: ${OUT_DIR}/${pdf_name}"
