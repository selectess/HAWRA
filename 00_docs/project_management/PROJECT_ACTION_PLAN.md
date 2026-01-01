# Plan d’Action Projet HAWRA

Version: 1.0
Date: 2025-12-12
Owner: Mehdi Wahbi (ORCID: 0009-0007-0110-9437)
DOI de référence: 10.5281/zenodo.17908061

## 1) Objectifs et Livrables
- O1. Réplicabilité immédiate (Docker): `docker run hawra` exécute "First Bloom" et produit `FIRST_BLOOM_REPORT.md`.
- O2. IDE Web Arbol (MVP): édition `H(q1)` dans le navigateur, compilation → visualisation Bloch/quantum_history.
- O3. Soumission scientifique: Nature/Science (Letter) avec figures et dépôt de code complet.
- O4. Préprint public: arXiv + Zenodo (liens actifs depuis le README).
- O5. Spécification BSIM v1.0 figée (schéma et compatibilité rétro maintenus).
- O6. Kit de diffusion: `HAWRA_Manifesto.md`, monographie HTML, Figure 1 (obsolescence cryogénique).

Livrables clés
- L1. `07_deployment/Dockerfile` (COMPLÉTÉ) + `docker-compose.yml` (optionnel) + README de lancement (à compléter).
- L2. `08_webapp` enrichi: éditeur Arbol, bouton Compiler/Simuler, rendu Bloch et logs.
- L3. `06_publication/HAWRA_Nature_Letter_Draft.md` finalisé + figures dans `06_publication/figures/`.
- L4. `README.md` (Badges, DOI, Quickstart 60s, contribution).
- L5. `00_docs/formalization/Project_Integrity_Nomenclature.md` stabilisé (versionné).

## 2) Ressources nécessaires
- Humaines
  - Responsable scientifique (PI) – arbitrage scientifique, soumission.
  - Ingénieur Dev/Compiler – Arbol/BSIM et simulateur.
  - Ingénieur DevOps – Docker, CI/CD, reproductibilité.
  - Designer/Comm. scientifique – figures, monographie, manifeste.
- Techniques
  - Python 3.12+, QuTiP, NumPy/Matplotlib; Docker 24+.
  - GitHub Actions (CI), GitHub Pages (docs) ou simple Flask local.
- Matérielles (phase numérique)
  - CPU standard; GPU non requis. (Phase wet-lab ultérieure non incluse ici.)

## 3) Calendrier et Jalons
- Semaine 1
  - J1: Dockerfile fonctionnel (L1) – First Bloom s’exécute en 60s.
  - J3: CI GitHub Actions (tests `pytest`, build image) – badge vert.
- Semaine 2
  - J6: IDE Web Arbol (MVP) – éditeur, compiler, afficher Bloch (L2).
  - J7: README Quickstart + Figures finalisées (L3/L4 partiel).
- Semaine 3
  - J12: Figée BSIM v1.0, nomenclature signée (L5).
  - J14: Dossier complet de soumission + préprint en ligne (L3/L4/L6).

Jalons d’acceptation
- M1: `docker run hawra` → "SUCCESS" dans `FIRST_BLOOM_REPORT.md`.
- M2: Compilation `.arbol` → `.bsim.json` via interface web et affichage de la sphère de Bloch.
- M3: CI passe (tests, lint), image Docker publiée.
- M4: Lettre Nature prête, figures référencées, DOI/ORCID cohérents.

## 4) Méthodologies et Outils
- Gestion de version: Git (branches `main` + `release/*`).
- Qualité: `pytest`, style PEP8, revues PR obligatoires.
- CI/CD: GitHub Actions (test → build → push image Docker Hub/GHCR).
- Documentation: Markdown + monographie HTML; issues/roadmap GitHub.
- Suivi technique: Kanban (To Do / In Progress / Done) et versionnage sémantique.

## 5) Contrôle Qualité et Validation
- Tests automatiques
  - `arbol/tests/*` (lexer/parser/compiler) – doivent passer à 100%.
  - `03_unified_simulator/tests/*` – stabilité BSIM et exactitude NOS/Simulateur.
- Contrôles de reproductibilité
  - Script non interactif: `scripts/run_first_bloom.py` – génère `FIRST_BLOOM_REPORT.md` et artefacts.
  - Image Docker reconstruite hebdomadairement via CI.
- Critères d’acceptation
  - Rapport "First Bloom" contient historique quantique, mesure et télémétrie biologique non vides.
  - Figure 1 comparant 300K vs 4K incluse et script reproductible.

## 6) Suivi et Reporting
- Rituels
  - Réunion hebdomadaire 30 min: avancées, blocages, décisions.
  - Point quotidien asynchrone (Slack/Issues) – micro-bilan et priorités.
- Indicateurs
  - Couverture tests, temps d’exécution simulation, taux de succès CI, téléchargements image Docker.
- Reporting
  - `00_docs/project_management/STATUS_REPORT.md` (hebdo): faits saillants, risques, prochaines actions.

## Annexes
- Quickstart (local)
  ```bash
  python -m venv .venv && source .venv/bin/activate
  pip install -r requirements.txt
  python scripts/run_first_bloom.py
  ```
- Quickstart (Docker)
  ```bash
  docker build -t hawra:latest .
  docker run --rm -it -v "$PWD":/workspace hawra:latest
  ```

---
Validation: ce plan doit être revu et approuvé en réunion de lancement avant exécution.

