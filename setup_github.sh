#!/bin/bash

# HAWRA Repository Setup Script
# Created by Move37 AI Team

echo "🌿 Initializing HAWRA Repository Structure..."

# Create core directories
mkdir -p src/bio-os
mkdir -p src/arbol-compiler
mkdir -p src/unified-sim
mkdir -p genomics/plasmids
mkdir -p genomics/sequences
mkdir -p docs/specs
mkdir -p docs/proofs
mkdir -p examples/arbol
mkdir -p examples/simulations
mkdir -p data/validation

# Create initial placeholder files
touch genomics/plasmids/.keep
touch genomics/sequences/.keep
touch data/validation/.keep

# Copy master files to root if they exist
[ -f GITHUB_README_MASTER.md ] && cp GITHUB_README_MASTER.md README.md
[ -f CITATION.cff ] && cp CITATION.cff CITATION.cff

# Create essential documentation files
cat <<EOF > CONTRIBUTING.md
# Contributing to HAWRA

We welcome contributions from the interdisciplinary community. 
Please follow the Code of Conduct and ensure all biological protocols 
comply with international safety standards.
EOF

cat <<EOF > SECURITY.md
# Security Policy

## Bio-Firewall & Digital Integrity
As HAWRA interacts with biological systems, any vulnerability that could 
compromise metabolic safety is treated as a critical security flaw.
EOF

echo "✅ HAWRA structure initialized successfully."
echo "🚀 Ready for git init and remote push."
