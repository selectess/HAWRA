# Format BSIM (Bio‑Simulation) — Arbol / HAWRA / PQPE

## Vue d’ensemble
- BSIM est le format de script produit par le compilateur Arbol pour piloter BioOS/HAWRA.
- Il encode des métadonnées et une séquence d’instructions normalisées.
- Par défaut, la sortie est un fichier JSON: `*.bsim.json`.
- L’extension `*.bsim` peut être utilisée (renommage), mais le contenu reste du JSON.

## Nom des fichiers
- Entrée: `*.arbol`
- Sortie: `*.bsim.json` (par défaut du compilateur).
- Option de déploiement: renommer en `*.bsim` si un moteur l’exige (même contenu).

## Structure du document
- `metadata`: informations sur la compilation (ex.: `source_arbol`, `version`).
- `instructions`: liste ordonnée d’instructions exécutables.

### Schéma des instructions
- `INITIALIZE`:
  - `config`:
    - `env.light_schedule`: liste de paires `(time, intensity)`.
    - `quantum`: paramètres (ex.: `p700_threshold`, `decoherence_rate`).
    - `bio`: paramètres (ex.: `p700_synthesis_rate`, `optimal_temp`).
- `QUANTUM_OP`:
  - `params.gate`: `H`, `X`, `Y`, `Z`, `CNOT` …
  - `params.qubits`: liste de noms de qubits.
- `MEASURE`:
  - `params.qubit`: nom du qubit.
  - `params.classical_bit`: nom du bit classique.
- `RUN_UNTIL`:
  - `params.time`: temps (secondes) d’exécution.
- `run`:
  - `circuit`: nom du circuit défini en Arbol.
  - `arguments`: dictionnaire de paramètres passés au circuit.
- `STIMULUS_APPLY` (générique):
  - `params.stimulus`: nom du stimulus.
  - `params.target`: cible (qubit ou entité).
  - `params.arguments`: paramètres d’application.

## Exemple minimal
```json
{
  "metadata": {"source_arbol": "compiled.arbol", "version": "0.1"},
  "instructions": [
    {
      "instruction_id": 0,
      "command": "INITIALIZE",
      "config": {
        "env": {"light_schedule": []},
        "quantum": {"p700_threshold": 0.8, "decoherence_rate": 0.05},
        "bio": {"p700_synthesis_rate": 1.0}
      }
    },
    {"instruction_id": 1, "command": "run", "circuit": "trainer", "arguments": {"lr": "0.1"}},
    {"instruction_id": 2, "command": "QUANTUM_OP", "params": {"gate": "H", "qubits": ["q1"]}},
    {"instruction_id": 3, "command": "MEASURE", "params": {"qubit": "q1", "classical_bit": "m"}}
  ]
}
```

## Génération (compilation)
- Depuis la racine du projet:
```
python3 -m 04_arbol.compiler.compiler 04_arbol/<script>.arbol
```
- Produit: `04_arbol/<script>.bsim.json`.

## Conventions et remarques
- Les valeurs d’arguments peuvent être des chaînes (ex.: `"1e-4"`) selon la grammaire actuelle.
- Le `light_schedule` est cumulatif si des stimuli de type `light` sont appliqués lors des runs.
- Les corps de circuits peuvent être déroulés à la compilation (configuration actuelle) pour produire `QUANTUM_OP`, `MEASURE`, `RUN_UNTIL`.

## Validation et analyse
- Le fichier `.bsim.json` est consommé par BioOS/HAWRA pour exécution et calcul des métriques.
- Visualisation possible via des scripts dédiés (ex.: comptage des instructions, trace du schedule lumineux).

## Évolutions prévues
- Option CLI pour choisir l’extension de sortie (`.bsim` vs `.bsim.json`).
- Paramètres continus et boucles contrôlées dans la grammaire Arbol.
- Schéma JSON formel (JSON Schema) pour validation automatique.
