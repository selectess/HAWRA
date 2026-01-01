# HAWRA Publication Methodology 2025

This document defines the standardized methodology for publishing research papers, creating complete dossiers, and managing the GitHub repository for the HAWRA project.

## 1. Research Paper Publication

### Standardized Compilation Process
1.  **Source Management**: All papers must be authored in LaTeX using the `templates/paper_template` structure.
2.  **Asset Linking**: Figures must be sourced from `06_presentation` or `results` and linked relative to the project root or copied via build script.
3.  **Compilation**:
    *   Primary: `pdflatex` or `lualatex` via CI/CD pipeline.
    *   Fallback: HTML generation with "Print-to-PDF" CSS for environments lacking TeX.
4.  **Output**: PDF (PDF/A compliant) and HTML (Web-ready).

### Quality Criteria
*   **Language**: Scientific English (US), checked via `textblob` or similar linter.
*   **References**: No missing citations (`[?]` or `??`). All bibtex entries must be valid.
*   **Figures**: Minimum 300 DPI for raster images. Vector graphics (SVG/PDF) preferred for plots.
*   **Reproducibility**: Every figure must link to the specific `simulation_run_id` and code version that generated it.

### Versioning
*   **Drafts**: `v0.x-alpha` (internal), `v0.x-beta` (collaborator review).
*   **Submission**: `v1.0-submission` (immutable tag).
*   **Revisions**: `v1.x-revision` (post-review).
*   **Final**: `v1.0-published` (DOI linked).

## 2. Complete Dossier Creation

### Structure
A "Complete Dossier" must follow this exact hierarchy:
```
HAWRA_Dossier_YYYY_MM_DD/
├── 00_Manifesto/          # Executive Summary & Vision
├── 01_Manuscript/         # The primary scientific paper (PDF + Source)
├── 02_Code/               # Snapshot of relevant source code (Arbol, BioOS)
├── 03_Data/               # Raw simulation data & GenBank files
├── 04_Multimedia/         # High-res figures & Videos/GIFs
└── 05_Supplementary/      # SOPs, Lab Protocols, Notebooks
```

### Mandatory Elements
*   `README.md`: Explaining the dossier contents.
*   `manifest.json`: List of all files with MD5 checksums.
*   `LICENSE`: Usage rights.
*   `requirements.txt` or `environment.yml`: For code reproducibility.

### Integrity Verification
*   Automated script `scripts/verify_dossier.py` must be run before distribution.
*   Checks: File existence, Checksum validity, PDF readability, Code syntax (linting).

## 3. GitHub Repository Management

### Architecture
Refactored structure for clarity and industry standards:
```
/
├── .github/               # CI/CD Workflows
├── src/                   # Source Code
│   ├── arbol/
│   ├── bioos/
│   └── simulator/
├── docs/                  # Documentation (Sphinx/MkDocs)
├── data/                  # Static Data (GenBank, Configs)
├── experiments/           # Jupyter Notebooks & One-off scripts
├── publications/          # LaTeX sources for papers
├── results/               # Simulation outputs (GitIgnored mostly)
├── scripts/               # Build & Maintenance utilities
└── tests/                 # Unit & Integration tests
```

### Conventions
*   **Naming**: `snake_case` for files/folders. `PascalCase` for Classes.
*   **Commits**: Conventional Commits (e.g., `feat: add new quantum gate`, `fix: correct plasmid sequence`).
*   **Branching**: `main` (stable), `develop` (integration), `feature/*` (dev).

### CI/CD Workflows
1.  **Test Suite**: Runs `pytest` on every push to `develop` and `main`.
2.  **Linting**: `flake8` or `black` for code quality.
3.  **Doc Build**: Generates static site from `docs/` to GitHub Pages.
4.  **Paper Compile**: Compiles LaTeX in `publications/` to PDF artifacts on release tags.

---

## Implementation Plan

1.  **Templates**: Create reusable templates for LaTeX papers and Dossier structures.
2.  **Scripts**: Develop `verify_publication.py` and `build_dossier.py`.
3.  **Migration**: Move existing files to the new `src/`, `experiments/`, `publications/` structure (Iterative process).
