#!/bin/bash
set -e
PROJECT_ROOT="/Users/mehdiwhb/Desktop/HAWRA"
BUNDLE_DIR="$PROJECT_ROOT/HAWRA_SUBMISSION_BUNDLE_2025_12_21"

echo "Creating bundle structure..."
mkdir -p "$BUNDLE_DIR/00_Manifesto"
mkdir -p "$BUNDLE_DIR/01_Manuscript"
mkdir -p "$BUNDLE_DIR/02_Code"
mkdir -p "$BUNDLE_DIR/03_Data"
mkdir -p "$BUNDLE_DIR/04_Multimedia"
mkdir -p "$BUNDLE_DIR/05_Supplementary"

echo "Copying Manifesto..."
cp "$PROJECT_ROOT/06_publication/HAWRA_Manifesto.md" "$BUNDLE_DIR/00_Manifesto/"
cp "$PROJECT_ROOT/06_publication/HAWRA_Nature_Submission.md" "$BUNDLE_DIR/00_Manifesto/"

echo "Copying Manuscript..."
cp "$PROJECT_ROOT/06_publication/latex/main.tex" "$BUNDLE_DIR/01_Manuscript/"
cp "$PROJECT_ROOT/06_publication/latex/references.bib" "$BUNDLE_DIR/01_Manuscript/"
cp -r "$PROJECT_ROOT/06_publication/latex/sections" "$BUNDLE_DIR/01_Manuscript/"
cp -r "$PROJECT_ROOT/06_publication/latex/figures" "$BUNDLE_DIR/01_Manuscript/"
cp "$PROJECT_ROOT/06_publication/latex/main.pdf" "$BUNDLE_DIR/01_Manuscript/HAWRA_Scientific_Paper_Wahbi_2025.pdf"

echo "Copying Code..."
cp -r "$PROJECT_ROOT/arbol" "$BUNDLE_DIR/02_Code/"
cp -r "$PROJECT_ROOT/bioos" "$BUNDLE_DIR/02_Code/"
cp -r "$PROJECT_ROOT/02_arbol_interface/jetson_client" "$BUNDLE_DIR/02_Code/"

echo "Copying Data..."
cp "$PROJECT_ROOT/03_quantum_simulation/results/phytoqmmml_convergence.json" "$BUNDLE_DIR/03_Data/" || true

echo "Copying Multimedia..."
cp -r "$PROJECT_ROOT/06_publication/latex/figures" "$BUNDLE_DIR/04_Multimedia/"

echo "Copying Supplementary..."
cp "$PROJECT_ROOT/00_docs/scientific/SOP_Implementation_Guide.md" "$BUNDLE_DIR/05_Supplementary/" || true
cp "$PROJECT_ROOT/00_docs/project_management/HAWRA_Scientific_FAQ_Defense.md" "$BUNDLE_DIR/05_Supplementary/" || true

echo "Creating README..."
cat <<EOF > "$BUNDLE_DIR/README.md"
# HAWRA Submission Bundle
Date: $(date +%Y-%m-%d)
Conceived & Engineered by: Mehdi Wahbi (Director, Move37 Initiative)
Implemented by: Move37 AI Team (Move37 Initiative)
Contact: m.wahbi.move37@atomicmail.io
DOI: 10.5281/zenodo.17908061

## Contents
- 00_Manifesto: Core vision and submission letters.
- 01_Manuscript: Full LaTeX source for the scientific paper and final PDF.
- 02_Code: Implementation of ARBOL, BioOS, and Jetson Interface.
- 03_Data: Simulation results and validation data.
- 04_Multimedia: High-resolution figures and diagrams.
- 05_Supplementary: SOPs and technical FAQs.
EOF

echo "Bundle created successfully."
