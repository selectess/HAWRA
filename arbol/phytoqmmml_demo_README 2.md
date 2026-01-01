# Démo PhytoQMML (Arbol / HAWRA / PQPE)

- Script: `phytoqmmml_demo.arbol`
- But: algorithme d'auto-apprentissage supervisé métabolique phyto-quantique
- Systèmes: HAWRA, PQPE, BioOS, langage Arbol
?
## Contenu

- `phytoqmmml_trainer`: préparation, corrélations, stimulus adaptatif lumineux, mesure (pseudo-label)
- `phytoqmmml_classifier`: évaluation de l'état via mesures quantiques
- `config`: sections bio/quantum/training
- `run`: exécution séquentielle (train puis classify)

## Compilation

```
python3 compiler/compiler.py phytoqmmml_demo.arbol
```

Sortie: `phytoqmmml_demo.bsim.json`

## Visualisation

Script de visualisation: `../03_unified_simulator/src/visualize_phytoqmmml.py`

Génère `../03_quantum_simulation/results/phytoqmmml_summary.png`

## Multi-runs et convergence

- Orchestration Arbol uniquement: `phytoqmmml_multi_runs.arbol`
- Compilation:

```
python3 -m 04_arbol.compiler.compiler 04_arbol/phytoqmmml_multi_runs.arbol
```

- Le `.bsim.json` résultant contient pour chaque run: `QUANTUM_OP`, `MEASURE`, `RUN_UNTIL` et un `light_schedule` accumulatif.
- L’analyse et la visualisation de convergence se font côté moteur de simulation (BioOS/HAWRA), en consommant le `.bsim.json`.
