# Cartographie Scientifique de l'Écosystème HAWRA

Cette cartographie répertorie l'intégralité des composants, documents et médias de la codebase HAWRA, servant de base structurelle pour la rédaction du papier scientifique de 137 pages.

## 1. Fondements Théoriques & Conceptuels
| Document | Rôle | Contenu Clé |
| :--- | :--- | :--- |
| `00_docs/concepts/HAWRA_Ecosystem_Overview.md` | Vue d'ensemble | Architecture ARBOL→BSIM→BioOS, Équation de Lindblad, Pipeline PQPE. |
| `00_docs/concepts/PQPE_Math_Bio_Engineering.md` | Ingénierie Mathématique | Modèles de coût hybrides, cinétique de Hill, thermodynamique du Silica Shield. |
| `00_docs/concepts/Hardware_API_Schema.md` | Interface Physique | Protocole de communication BioOS ↔ Jetson (REST API). |

## 2. Génomique & Ingénierie Génétique (01_genomics)
| Document / Fichier | Rôle | Contenu Clé |
| :--- | :--- | :--- |
| `01_genomics/plasmids/validated/HAWRA_FINAL_VALIDATED.gb` | Blueprint Génétique | Séquence annotée du plasmide (CaMV 35S, Lsi1, psaA, LUC). |
| `01_genomics/raw_sequences/CDS/` | Librairie de Gènes | Séquences FASTA/GB pour chaque composant (psaA, LUC, NOS, etc.). |
| `01_genomics/genome_analysis_scripts/visualize_plasmid.py` | Visualisation | Script de génération des cartes circulaires du plasmide. |
| **Artifacts Visuels** | | |
| `05_data/results/genetic_visualization/hawra_plasmid_map.png` | Carte du Plasmide | Visualisation 2D des promoteurs et CDS. |
| `05_data/results/plasmid_pqpe_3d.gif` | Modélisation 3D | Rendu dynamique de la structure tertiaire de l'ADN. |

## 3. Langage de Programmation & Compilation (arbol)
| Document / Fichier | Rôle | Contenu Clé |
| :--- | :--- | :--- |
| `arbol/compiler/compiler.py` | Moteur de Compilation | Parser Lark v0.3 traduisant l'Arbol en instructions BSIM JSON. |
| `arbol/compiler/arbol_ast.py` | Structure de Données | Définition de l'Arbre de Syntaxe Abstraite (AST) pour les opérations bio-quantiques. |
| `arbol/test.arbol` | Benchmarking | Script de test couvrant stimuli, qubits logiques et mesures. |

## 4. Simulation Multiphysique & Unified Simulator (03_unified_simulator)
| Document / Fichier | Rôle | Contenu Clé |
| :--- | :--- | :--- |
| `03_unified_simulator/src/hawra_simulator/simulator.py` | Moteur d'Exécution | Intégration des moteurs Quantique, Biologique et Environnemental. |
| `03_unified_simulator/src/hawra_simulator/engines/quantum_core.py` | Noyau Quantique | Simulation des portes (H, X, CNOT) via QuTiP et Lindblad. |
| `03_unified_simulator/src/hawra_simulator/engines/biological_system.py` | Noyau Biologique | Dynamique de la chlorophylle P700 et régulation GRN. |
| **Artifacts Visuels** | | |
| `bioos/simulations/multiphysics_simulator/p700_bloch_sphere.gif` | État Quantique | Évolution du qubit sur la sphère de Bloch. |
| `05_data/pqpe_simulation/pqpe_dynamics_real_units.png` | Cinétique Bio | Courbes de concentration et facteurs de cohérence. |

## 5. Système d'Exploitation Biologique (bioos)
| Document / Fichier | Rôle | Contenu Clé |
| :--- | :--- | :--- |
| `bioos/core/bio_os.py` | Noyau de l'OS | Gestionnaire de ressources, orchestration hardware et logging. |
| `bioos/core/isolation_control.py` | Contrôle d'Isolation | Algorithmes de stabilisation thermique et électromagnétique. |
| `bioos/quantum_interface/readout_improved.py` | Acquisition de Données | Traitement du signal de sortie (Luciférase / Électrodes). |

## 6. Interface Cyber-Physique & Client (02_arbol_interface)
| Document / Fichier | Rôle | Contenu Clé |
| :--- | :--- | :--- |
| `02_arbol_interface/jetson_client/client.py` | Client Matériel | Serveur Flask pour le pilotage des actionneurs sur NVIDIA Jetson. |
| `08_webapp/app.py` | Monitoring Web | Dashboard temps-réel et visualisation des résultats de simulation. |

## 7. Résultats & Preuves de Concept
| Document / Fichier | Rôle | Contenu Clé |
| :--- | :--- | :--- |
| `05_data/results/simulation_log.json` | Données Brutes | Logs complets d'une exécution PQPE. |
| `results/hawra_poc_execution_v2.gif` | Preuve d'Exécution | Démonstration visuelle du cycle complet ARBOL → Hardware. |
| `06_publication/figures/figure1_metabiotic_advantage.png` | Analyse Scientifique | Graphique comparatif de l'avantage PQPE vs Silicium. |

## 8. Workflows & Pipelines
1.  **Workflow de Conception Génétique:** `raw_sequences` → `CDS selection` → `validated GB` → `visualize_plasmid`.
2.  **Workflow de Compilation:** `.arbol` → `Lark Parser` → `.bsim.json`.
3.  **Workflow d'Exécution:** `BioOS` → `Simulator` (Digital Twin) ↔ `Jetson Client` (Physical Reality).
4.  **Workflow de Visualisation:** `simulation_log.json` → `WebApp` → `3D Renders`.
