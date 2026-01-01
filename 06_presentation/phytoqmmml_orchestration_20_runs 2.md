# PhytoQMML Orchestration (20 runs) — HAWRA/PQPE/Arbol

## Objectif
- Démontrer une orchestration multi-runs entièrement en Arbol (20 runs) pour PhytoQMML.
- Dérouler les corps des circuits à la compilation (`QUANTUM_OP`, `MEASURE`, `RUN_UNTIL`).
- Produire un `light_schedule` cumulatif pour refléter les stimuli lumineux dans le temps.

## Scripts et compilation
- Script Arbol: `04_arbol/phytoqmmml_multi_runs.arbol`
- Compiler:
```
python3 -m 04_arbol.compiler.compiler 04_arbol/phytoqmmml_multi_runs.arbol
```
- Sortie:
  - `04_arbol/phytoqmmml_multi_runs.bsim.json`
  - Contient pour chaque run: opérations quantiques, mesures, et `RUN_UNTIL`.

## Orchestration
- 20 runs d’entraînement (`phytoqmmml_trainer`) alternés avec un circuit `phytoqmmml_monitor` (multiples mesures) pour affiner l’état.
- `light_schedule` cumulatif (exemple): `0→10→20→…→200 s` avec alternance intensité/arrêt.

## Méthodologie et analyse
- Le `.bsim.json` est consommé par BioOS/HAWRA pour calculer les métriques (fidelité, perte, etc.).
- Recommandé:
  - Visualiser les comptages (`QUANTUM_OP`, `MEASURE`, `RUN_UNTIL`).
  - Tracer l’évolution de métriques par run (produites côté BioOS/HAWRA).

### Sorties métriques (auto-générées par BioOS)
- `05_data/results/simulation_log.json`
- `05_data/results/bsim_metrics.json` (par run: complexité, fidélité proxy, perte proxy)
- `05_data/results/bsim_convergence.png` (courbes fidélité/perte)

## Extraits (référence Arbol)
- Stimulus
```
stimulus adaptive_light(intensity: float, wavelength: float, duration: float, type: light) { }
```
- Entraînement
```
run phytoqmmml_trainer(learning_rate: 0.15, coherence: 1e-4, noise_level: 0.05);
```
- Monitoring
```
run phytoqmmml_monitor(samples: 4);
```

## Points clés
- Transparence et auditabilité: orchestration en Arbol, inspection via `.bsim.json`.
- Compatibilité PQPE/HAWRA: déroulement des circuits et schedule lumineux intégrés.
- Base solide pour la convergence supervisée: horizon étendu (20 runs) et mesures denses.

## Prochaines étapes
- Étendre le monitor (plus de qubits/mesures) pour une granularité accrue.
- Ajouter des runs de classification intermédiaires (`phytoqmmml_classifier`) pour suivre les états.
- Préparer figures de convergence basées sur les métriques BioOS/HAWRA.
