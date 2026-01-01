t# PQPE – Documentation Mathématique, Phyto‑Ingénierie Génétique et Computationnel

## Résumé

HAWRA met en œuvre une PhytoQuantum Processing Entity (PQPE): une plante synthétique programmable qui exécute des opérations quantiques via des stimuli physiques (lumière/EM) traduits par ARBOL et interprétés par un BioOS vivant. Ce document formalise le cadre mathématique, la phyto‑ingénierie génétique, l’architecture computationnelle et les critères de validation numérique.

## Objectifs et Critères

- Objectif: décrire et valider un système PQPE end‑to‑end (plasmide → ARBOL/BSIM → simulation BioOS+quantique → métriques/figures → bundles).
- Critères de réussite:
  - Génération d’un GenBank “zéro placeholder” avec promoteurs/terminators réels (CaMV 35S/NOS).
  - Compilation et exécution de scénarios ARBOL en BSIM sans erreurs.
  - Observables quantiques ⟨Z⟩ enregistrées et tracées.
  - Métriques consolidées (cohérence moyenne, P700, lumière, score énergétique).
  - Figures comparatives et spectrales, bundles préprint/archives.

## Modèle Mathématique

- État du qubit biologique (P700): espace de Hilbert à deux niveaux, `|ψ⟩ = α|0⟩ + β|1⟩` avec `|α|²+|β|²=1`.
- Évolution (Lindblad open‑system): `dρ/dt = -i [H(t), ρ] + Σ_k γ_k (L_k ρ L_k^† - ½{L_k^† L_k, ρ})`.
  - Hamiltonien de contrôle: `H(t) = H_0 + Ω(t) σ_x + Δ(t) σ_z`.
  - Déphasage et relaxation: `L_T2* = √(1/T2*)·σ_z`, `L_T1 = √(1/T1)·σ⁻`.
- Couplage spectral (lumière): `κ(λ) = exp(−(λ−λ_peak)²/(2σ²))`, intensité effective: `I_eff(t) = I(t)·κ(λ)` → `Ω(t) = α·I_eff(t)`.
- Mesure ensemble: `⟨Z⟩_ensemble(t) = (1/N) Σ_i Tr(ρ_i(t) σ_z)`; lecture bioluminescente par ratio vert/rouge `R = I_verte/(I_verte+I_rouge)`.
- Score énergétique: `score = coherence_avg + p700_avg − 0.01 · light_avg` (agrégation numérique).

## Phyto‑Ingénierie Génétique

- Cassette computationnelle modulaire par gène: Promoteur (CaMV 35S) → Insulator (RiboJ) → RBS (B0034) → CDS (psaA/CRY2/Luc/Lsi1/HSP70/PEPC) → Terminator (NOS).
- Génération GenBank:
  - Script: `scripts/create_genbank.py` (préférence FASTA 35S/NOS; labels `P35S_`/`NOS_`).
  - Sortie: `05_data/results/hawra_plasmid.gb`.
- Rôles: psaA (qubit P700), CRY2 (portes), Luc verte/rouge (lecture), Lsi1 (cage silice), HSP70 (thermo‑protection), PEPC (énergie CAM).

### BioOS Map — Correspondance Gène → Fonction OS

| Gène | Accession (ex.) | Rôle dans BioOS | Mécanisme / Notes |
|---|---:|---|---|
| `psaA` | `NC_031161.1` | Qubit cœur (registre) | Centre réactionnel P700; état quantique (|0>,|1>) lié à excitons photosynthétiques; cible de lecture via luciferase couplée |
| `CRY2` | `AF156319.1` | Porte quantique / détecteur d’instruction | Photoreceptor sensible au bleu et EM; conformation controlée par stimulus → interaction avec P700 (appliquer rotations) |
| `CRY1` | `NM_119628` | Porte stable / modulation | Photoreceptor secondaire pour biais probabiliste et gating temporel |
| `Luc` (verte) | `U47132.1` | Readout (|0>) | Riboswitch ATP → expression Luc verte (560 nm) signale état stable; utilisé pour lecture optique |
| `Luc` (rouge, S284T) | `U47132.1 (S284T)` | Readout (|1>) | Riboswitch ROS → expression Luc rouge (615 nm) signale état instable |
| `Lsi1` / `SIT1` | `KT159333.1` | Isolation quantique / cage de silice | Transporteur silicium → déposition contrôlée de silice autour chloroplastes; augmente T2 (+66.67%) |
| `HSP70` ×2 | `M60185` | Gestionnaire d’erreurs thermique / maintenance | Chaperon protéique, protège l’intégrité structurelle sous stress thermique; active réparation cellulaire |
| `PEPC` | `X64138` | Gestion énergétique (CAM) | Module métabolique CAM pour approvisionnement ATP nocturne et découplage jour/nuit |
| `dCas9 + gRNA` | (tooling) | Mémoire épigénétique / apprentissage | Méthylation/épimodification ciblée des promoteurs; modulable par stimuli répétés pour adaptation transgénérationnelle |

Cette table explicite la correspondance entre les composants génétiques présents dans la cassette HAWRA et leurs fonctions d’OS biologique. Elle sert de référence pour le design expérimental, les tests de validation et la traçabilité des comportements observés en simulation et en laboratoire.

## Architecture Computationnelle

- Langage ARBOL:
  - Grammaire: `arbol/grammar/arbol_unified.ebnf`.
  - Compilateur: `arbol/compiler/compiler.py` (instructions BSIM; stimuli `visit_StimulusApplication`).
  - Programmes démo: `arbol/*.arbol`.
- BioOS (simulation vivante):
  - Environnement (température, lumière, wavelength): `03_unified_simulator/src/hawra_simulator/engines/environment.py:15-27,57-75,76-83`.
  - Biologie (GRN+Hill+spectral): `03_unified_simulator/src/hawra_simulator/engines/biological_system.py:87-110`.
  - Stimuli et planning: `03_unified_simulator/src/hawra_simulator/simulator.py:203-282`.
  - Observables ⟨Z⟩ (QuTiP): `03_unified_simulator/src/hawra_simulator/simulator.py:321-334`.
- Exploitation:
  - Agrégateur: `03_unified_simulator/aggregate_outputs.py` (CSV/JSON, score, parité schedule/instructions).
  - Figures: `03_unified_simulator/src/plot_comparisons.py` (comparatifs, comptages ⟨Z⟩, `gene_spectral_response`, sweep spectral CSV/PNG).
  - CLI one‑shot: `scripts/one_shot_pipeline.py` (compile→simulate→aggregate→figures→bundles).

## Protocoles et Validation

- Protocole expérimental: `00_docs/concepts/experimental_validation_plan.md` (+ `experimental_validation_plan.rtf`).
- Étapes numériques:
  - Compiler ARBOL: `./.venv/bin/python -m arbol.compiler.compiler arbol/pqpe_simulation.arbol`.
  - Exécuter BSIM: `./.venv/bin/python dist/preprint_submission/scripts/run.py --bsim 03_unified_simulator/pulses_40_1_5.bsim.json --output 03_unified_simulator/output_pulses_40_1_5.json`.
  - Agréger métriques: `./.venv/bin/python 03_unified_simulator/aggregate_outputs.py --pattern '03_unified_simulator/output_*.json' --format csv --out 03_unified_simulator/results/consolidated_metrics.csv --score --sort_by score --desc`.
  - Générer figures: `./.venv/bin/python 03_unified_simulator/src/plot_comparisons.py`.
  - Pipeline one‑shot: `./.venv/bin/python scripts/one_shot_pipeline.py --bsim 03_unified_simulator/pulses_40_1_5.bsim.json --out 03_unified_simulator/results/one_shot_output.json`.

## Résultats et Indicateurs

- Configuration gagnante: `pulses_40_1_5` (cohérence ~0.495, score max), efficacité énergétique supérieure aux références.
- Réponse spectrale: gA(650 nm) > gB(450 nm) sous λ=650; sweep `spectral_sweep.csv/png`.
- Observables ⟨Z⟩: `phytoqmmml_quantum_counts_observables.png`.
- Bundles: préprint et archive complète (`dist/preprint_submission/…zip`, `dist/full_archive/…zip`).

## Risques et Pistes d’Amélioration

- Séquences exactes 35S/NOS: vérifier/replacer FASTA pour GenBank “zéro placeholder”.
- Modélisation EM/CRY2: affiner couplage et paramètres avec mesures in vivo.
- Intégration Jetson I/O: définir JSON de pilotage (LED/EM) et watchdog sécurité (T/ROS).
- Tests paramétriques stimuli combinés et extension des observables.

## Conclusion

Le cadre PQPE de HAWRA est formalisé et validé numériquement: math/physique, génétique, computation, simulation vivante et observables quantiques. Les artefacts garantissent la reproductibilité et préparent la phase expérimentale. Les prochaines intégrations (FASTA 35S/NOS exacts, I/O Jetson) mèneront à une mise en culture et à la validation in vivo.