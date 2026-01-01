# Nomenclature Révisée et Normalisation (HAWRA/PQPE/BioOS/Arbol)

Ce document révisé définit la nomenclature, les conventions et les schémas de données standard pour l’ensemble du projet HAWRA. Il vise à:
- Normaliser les noms et structures (DSL Arbol, BSIM, simulateur, BioOS, résultats).
- Garantir la compatibilité entre code, scripts et documents.
- Servir de référence unique pour le développement et l’expérimental.

## Portée et Références
- Périmètre: `Arbol` (langage/compilation), `BSIM` (assembly JSON), `03_unified_simulator` (exécution unifiée), `BioOS` (multiphyse), données génomiques et webapp.
- Références: `PhytoQMML_formalization.md`, `HAWRA-PQPE_formal_model.md`, `Hardware_API_Schema.md`, exemples `.arbol` et `.bsim.json`.

## Conventions Générales
- Casse et style:
  - Clés JSON et champs de config: `snake_case`.
  - Noms de portes quantiques: `UPPERCASE` (`H`, `X`, `Y`, `Z`, `CNOT`, `CX`).
  - Identifiants de gènes/qubits/bits: `lower_snake_case` ou `q1`, `c1` pour simplicité.
- Types et unités:
  - `duration`: secondes (float). Chaînes suffixées acceptées par le compilateur (`ms`, `s`, `min`, `h`) et normalisées en secondes.
  - `wavelength`: nanomètres (float).
  - `intensity`: sans unité (float), interprétée comme intensité relative (0..∞). Les planchers/plafonds sont définis par scénario.
  - `temperature`: °C (float).
- Fichiers:
  - Sources DSL: `.arbol`.
  - Assembly: `.bsim.json`.
  - Résultats: `.json/.csv/.png/.gif`.
  - Génomique: `.gb`, `.fasta`.

## Concepts Paradigmatiques (Metabiotic Terminology)
- **Metabiotic Architecture** : Paradigme informatique où le matériel de calcul est un organisme vivant génétiquement modifié, par opposition au hardware inerte (silicium/cryogénie).
- **PhytoQMMML** (Phyto-Quantum Math-Metabolic Machine Learning) : Cadre d'apprentissage automatique supervisé utilisant les flux métaboliques comme vecteurs d'optimisation et la cohérence quantique comme fonction de coût.
- **Bio-Qubit** : Unité d'information quantique portée par une structure biologique (ex: P700 dans le Photosystème I).
- **Phyto-Programming** : Acte de coder des fonctions logiques et quantiques s'exécutant sur le métabolisme végétal via le DSL ARBOL.
- **BQG-IR** (Bio-Quantum Graph Intermediate Representation) : Représentation intermédiaire du compilateur ARBOL permettant d'optimiser les pulses photoniques en fonction de l'état métabolique.
- **GRAPE Optimization** : Algorithme utilisé par le compilateur pour sculpter l'enveloppe des pulses lasers afin de maximiser la fidélité des portes.
- **Bio-SGD** : Algorithme de descente de gradient stochastique où les poids synaptiques sont représentés par des niveaux de méthylation de l'ADN.
- **Silica Shield** (Lsi1) : Mécanisme de protection contre la décohérence via l'expression du gène `Lsi1`, favorisant la biominéralisation d'une structure de silice autour des thylakoïdes.
- **Wetware-Reliant** : Qualifie une architecture qui ne peut fonctionner sans substrat biologique actif.

## Structure du Dossier de Soumission (Submission Bundle)
Conformément à `METHODOLOGY_2025.md`, tout dépôt officiel doit suivre cette structure :
- `00_Manifesto/` : Résumé exécutif et vision (Mehdi Wahbi, Directeur Move37 Initiative).
- `01_Manuscript/` : Manuscrit scientifique (LaTeX + PDF/A).
- `02_Code/` : Snapshots de `arbol/` et `bioos/`.
- `03_Data/` : Données de simulation (`.json`, `.csv`) et fichiers génomiques (`.gb`).
- `04_Multimedia/` : Figures haute résolution et schémas TikZ/Mermaid.
- `05_Supplementary/` : Protocoles SOP (SOP-01 à SOP-06).

## Auteur et Crédits
- **Concepteur & Architecte Principal** : Mehdi Wahbi (Directeur de Move37 Initiative).
- **Implémentation & Validation** : Move37 AI Team (Move37 Initiative).
- **ORCID** : 0009-0007-0110-9437.
- **Contact** : m.wahbi.move37@atomicmail.io
- **Citation** : *"On ne rêve plus, on compile."*

## Schéma BSIM Standard (instructions)
- `INITIALIZE`:
  - Champ: `config` (objet) avec sous‑sections:
    - `dt` (float), `max_time` (float optionnel).
    - `env`: voir section Environnement.
    - `bio`: voir section Biologie.
    - `quantum`: voir section Quantique.
- `RUN_UNTIL`:
  - `params.time` (float, en secondes).
- `STIMULUS_APPLY`:
  - `params.stimulus` (string): `light`, `adaptive_light`, `heat`, `chemical`, `inhibitor`…
  - `params.target` (string|null) selon le stimulus.
  - `params.arguments` (objet): clés normalisées; ex. `intensity`, `duration`, `wavelength`, `temperature`, `effect`, `magnitude`.
- `QUANTUM_OP`:
  - `params.gate` (string): `H|X|Y|Z|CNOT|CX`.
  - `params.qubits` (array): liste d’indices (int) ou noms (string) des qubits.
- `MEASURE`:
  - `params.qubit` (int|string).
  - `params.classical_bit` (int|string|null).

### Normalisation et compatibilité
- Synonymes acceptés par le chargeur (normalisés à l’exécution):
  - `op_code`/`opcode` → `command`.
  - `target`/`targets` → `qubits` (pour `QUANTUM_OP`), sinon `target` reste tel quel pour `STIMULUS_APPLY`.
  - `cbit`/`classical_bit` → `classical_bit`.
  - `params` absent → `{}`.
- Les portes multi‑qubits supportées explicitement: `CNOT` (ou `CX`). Autres portes 2‑qubits requièrent une extension.

## Environnement (config.env)
- Champs standard:
  - `temperature` (float, °C).
  - `light_intensity` (float, par défaut 0.0).
  - `light_wavelength` (float|null).
  - `light_schedule` (array d’événements):
    - Événement: `{ time: float, intensity: float, wavelength?: float }`.
  - Comportements:
    - Ordonnancement trié par `time` et compactage (fusion des plateaux) recommandé.
    - Curseur interne `_cursor` géré par l’environnement.

## Biologie (config.bio)
- Gènes:
  - `genes`: liste d’objets gène:
    - `{ id: string, basal_rate: float, degradation_rate: float, light_sensitivity?: float, peak_wavelength?: float, peak_sigma?: float, initial_expression?: float }`.
- Réseau de régulation (GRN):
  - `grn`: dictionnaire `{ target_gene: [reg_info, ...] }`.
  - `reg_info`:
    - `{ id: string, type: activator|repressor, weight: float, hill_coefficient: float, half_max_concentration: float }`.
  - Ancien format (dict `{regulator: weight}`) converti en interne vers la forme normalisée.
- Paramètres globaux:
  - `synthesis_rate` (float), `initial_p700` (float), `optimal_temp` (float), `temp_sigma` (float).
- Configuration Entraînement (PhytoQMMML):
  - `config.training`:
    - `learning_rate` (float): Taux d'apprentissage pour l'adaptation métabolique.
    - `target_efficiency` (float): Seuil de rendement quantique visé.

## Quantique (config.quantum)
- Paramètres physiques:
  - `t1` (float, s), `t2` (float, s), `silica_protection_factor` (float).
- États de base:
  - `psi0` implicite (|0⟩), configurable ultérieurement si nécessaire.
- Simulation:
  - Hamiloniens de contrôle via `control_pulse { amplitude }`.
  - Opérateurs de bruit: relaxation (T1), déphasage (T2).

## Nomenclature Arbol (DSL)
- Mots clés: `circuit`, `gate`, `stimulus`, `apply`, `measure`, `config`, `run`, `LOGICAL_QUBIT`, `CLASSICAL_BIT`.
- Portes natives: `H`, `X`, `Y`, `Z`, `CNOT` (`CX` en alias pour compilation).
- Stimulus `light`/`adaptive_light`: arguments normalisés (`intensity`, `duration`, `wavelength`).
- Mesure: `measure q [ON|->|TO] c;` (bit classique optionnel).
- Configuration: `config { bio { ... }, quantum { ... }, env { ... } }`.

## Conventions de Nommage Fichiers/Dossiers
- Dossiers fonctionnels: 
  - Langage/compilation: `arbol/`.
  - Simulateur unifié: `03_unified_simulator/` (et alias `bioos/simulations/multiphysics_simulator/`).
  - BioOS: `bioos/`.
  - Simulation Quantique Spécifique: `03_quantum_simulation/`.
  - Données: `05_data/`, `results/`, `01_genomics/`.
  - Publication & Figures: `06_publication/`.
  - Documentation Scientifique: `00_docs/scientific/`.
- Fichiers:
  - Scripts Python: `lower_snake_case.py`.
  - Docs: `*.md`, `*.tex`.
  - Sorties: `*.json`, `*.png`, `*.gif`, `*.csv`.

## Gestion des Données et Artefacts (Granularité)
- **Logs de Simulation :** Stockés dans `results/` ou `bioos/simulations/logs/`.
  - Format: `simulation_name_timestamp.json` ou `simulation.log`.
- **Figures Générées :**
  - Plots statiques: `*.png` (ex: `coherence_plot.png`, `gene_regulation_p700.png`).
  - Animations: `*.gif` (ex: `bloch_sphere_decoherence.gif`).
- **Checkpoints de Données :**
  - `validation_results.json`: Résultats agrégés des tests de validation.
  - `bsim_metrics.json`: Métriques de performance du compilateur/simulateur.

## Compatibilité et Migration
- Chargeur BSIM assure la compatibilité avec anciens schémas (`opcode`, GRN dict, cbit list…).
- Pour nouveaux scripts, utiliser les clés normalisées de ce document.
- Les tests assurent la parité `light_schedule` vs `STIMULUS_APPLY` (`03_unified_simulator/tests/test_simulator.py`).

## Exemples Minimaux
### Stimulus lumineux adaptatif
```
stimulus adaptive_light(intensity: float, wavelength: float, duration: float, type: light) {}
apply adaptive_light(intensity: 0.8, wavelength: 680.0, duration: 10.0) on q1;
```
### Circuit et exécution
```
circuit trainer(learning_rate: float, coherence: float) {
  q1: qubit; q2: qubit; m: bit;
  H(q1); CNOT(q1, q2);
  m = measure q2;
}
run trainer(learning_rate: 0.15, coherence: 1e-4);
```

## Checklist de Validation
- Clés JSON en `snake_case` et unités conformes.
- `INITIALIZE` contient `env`, `bio`, `quantum`, `dt`.
- Schedules triés et compactés, événements `{time,intensity,wavelength?}`.
- GRN au format liste de régulateurs avec paramètres Hill.
- Portes multi‑qubits limitées à `CNOT` tant que non étendues.
- Tests de parité et de décohérence passent (`pytest`).

## Conclusion
Cette révision unifie la terminologie et la structure des données à travers Arbol, BSIM, le simulateur et BioOS. Elle facilite la programmation du « qubit vivant », l’intégration expérimentale et la maintenance logicielle, tout en préservant la compatibilité avec les artefacts existants.
