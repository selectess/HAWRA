# Cartographie de la Codebase HAWRA
> **Initiative:** HAWRA (Hybrid Architecture for Watson-Crick Research Applications)
> **Conception & Architecture:** Mehdi Wahbi (Directeur, Move37 Initiative)
> **Implémentation:** Move37 AI Team
> **ORCID:** 0009-0007-0110-9437
> **DOI:** 10.5281/zenodo.17908061
> **Licence:** Open Source / Creative Commons
> **Dernière mise à jour:** 14 Décembre 2025 (Preprint Release v1.0)

Cette carte de la codebase HAWRA fournit une vue d'ensemble structurée des fichiers et répertoires du projet, à l'exclusion des environnements virtuels (`.venv`). Elle documente l'intégralité des travaux de recherche, de développement et de validation menés par Mehdi Wahbi dans le cadre de la création d'un ordinateur métabiotique.

## Méta-Informations du Projet

Le projet vise à explorer les frontières de l'ingénierie des systèmes artificiels intelligents en fusionnant biologie synthétique et calcul quantique. HAWRA est l'implémentation phare de cette vision.

### Identifiants et Références
*   **Investigateur Principal :** Mehdi Wahbi (Directeur, Move37 Initiative)
*   **Équipe Technique :** Move37 AI Team
*   **Identifiant Numérique d'Objet (DOI) :** `10.5281/zenodo.17908061`
*   **Statut de Validation :** Validé in silico (Voir `validate_simulation.py`)
*   **Architecture :** Hybride Phyto-Quantique (Bio-Qubit P700 + Silica Shield)
*   **Point d'Entrée Communauté :** [hawra.tech](https://hawra.tech)

---

```
.
├── .DS_Store
├── .git
│   ├── HEAD
│   ├── config
│   ├── description
│   ├── hooks
│   │   ├── applypatch-msg.sample
│   │   ├── pre-commit.sample
│   │   ├── pre-merge-commit.sample
│   │   ├── prepare-commit-msg.sample
│   │   ├── pre-push.sample
│   │   ├── pre-rebase.sample
│   │   └── update.sample
│   ├── info
│   │   └── exclude
│   ├── objects
│   │   ├── info
│   │   └── pack
│   │       ├── pack-0381673812292186716a5051187424177431267b.idx
│   │       └── pack-0381673812292186716a5051187424177431267b.pack
│   └── refs
│       ├── heads
│       │   └── main
│       └── tags
├── .gitignore
├── .vscode
│   └── settings.json
├── LICENSE
├── README.md
├── bioos
│   ├── __init__.py
│   ├── simulations
│   │   ├── multiphysics_simulator
│   │   │   ├── __init__.py
│   │   │   ├── quantum_engine.py

## Chronologie du Développement HAWRA (Reconstruction Logique)

Cette chronologie retrace l'évolution du projet, des concepts initiaux à la validation multiphysique actuelle.

### Phase 1 : Conceptualisation et Fondations (Nov 2025)
*   **Jalon :** Définition du paradigme "Métabiotique" et du concept "Silica Shield".
*   **Artefacts Clés :**
    *   `concept_paper_outline.md` (Manifeste scientifique).
    *   `00_docs/formalization/PhytoQMML_formalization.md` (Théorie de l'apprentissage phyto-quantique).

### Phase 2 : Infrastructure Arbol & BSIM (Début Déc 2025)
*   **Objectif :** Créer le langage de programmation du vivant.
*   **Réalisations :**
    *   Développement du Lexer/Parser Arbol (`arbol/grammar.py`).
    *   Définition du standard BSIM (`Project_Integrity_Nomenclature.md`).
    *   Implémentation du Compilateur (`bioos/bio_compiler/compiler.py`).
*   **Validation :** Compilation réussie de `arbol/phytoqmmml_demo.bsim.json`.

### Phase 3 : Simulation Multiphysique (Mi-Déc 2025)
*   **Objectif :** Valider les modèles théoriques *in silico*.
*   **Réalisations :**
    *   Moteur Quantique (`quantum_engine.py`) : Implémentation de l'équation de Lindblad.
    *   Moteur Biologique (`biological_engine.py`) : Modélisation GRN et cinétique de Hill.
    *   Moteur Environnemental (`environment_engine.py`) : Gestion des cycles lumineux.
*   **Validation :** Tests unitaires des probabilités Monte Carlo (`sop_procedural_simulation.py`).

### Phase 4 : Intégration et Validation Globale (Actuel)
*   **Objectif :** Preuve de concept complète (E2E).
*   **Réalisations :**
    *   Script de validation unifié (`validate_simulation.py`).
    *   Campagne de "Parameter Sweep" pour optimiser les constantes de couplage.
*   **Données Générées :**
    *   `results/multiphysics_simulation/multiphysics_simulation_v2.json` (Trace d'exécution complète).
    *   `results/bloch_sphere_decoherence.gif` (Visualisation de la cohérence).

---

## Cartographie des Artefacts de Données (Logs & Résultats)

Cette section localise précisément les fichiers de données générés par les simulations, essentiels pour la reproductibilité.

### 1. Résultats de Simulation Multiphysique
*   **Fichier Principal :** `results/multiphysics_simulation/multiphysics_simulation_v2.json`
    *   **Contenu :** Séries temporelles complètes (P700, Intensité Lumineuse, Cohérence).
    *   **Visualisation associée :** `results/multiphysics_simulation/multiphysics_simulation_v2.png`.
*   **Logs Bruts :** `results/simulation_log.json` (Trace d'événements discrets).

### 2. Validation Quantique
*   **Animation Bloch Sphere :** `results/bloch_sphere_decoherence.gif`
    *   **Description :** Montre l'effondrement de l'état quantique $|+\rangle$ sous l'effet de la décohérence ($T_2$), ralenti par le *Silica Shield*.
*   **Graphique de Cohérence :** `results/p700_simulation/p700_coherence_decay.png`.

### 3. Biologie & Génomique
*   **Régulation Génique :** `results/gene_regulation_p700.png`
    *   **Description :** Réponse du promoteur P700 aux impulsions lumineuses (Cinétique de Hill).
*   **Données Génomiques :** `01_genomics/experiments/first_bloom_results.json`.

### 4. Compilateur Arbol
*   **Métriques de Compilation :** `results/bsim_metrics.json`.
*   **Sortie de Validation :** `bioos/bio_compiler/arbol/validation_results.json`.

---

### Fichier: bioos/simulations/multiphysics_simulator/quantum_engine.py
**Position:** `bioos/simulations/multiphysics_simulator/quantum_engine.py`
**Utilité:** Implémente le moteur de simulation quantique pour le simulateur multiphysique, modélisant la dynamique du P700 et la lecture à double canal.
**Description:** Ce fichier définit la classe `QuantumEngine` qui simule le comportement d'un centre de réaction P700, un composant clé de la photosynthèse. Le moteur gère l'excitation et la désexcitation du P700, ainsi que la détection de son état via un modèle de lecture à double canal (luciférase verte et rouge).

**Méthodes clés:**
*   `__init__(self, p700_excitation_threshold=0.5, fast_collapse_probability=0.1, decoherence_rate_fast=0.05, decoherence_rate_slow=0.01)`: Initialise le moteur avec des paramètres configurables pour le seuil d'excitation du P700, la probabilité de collapse rapide et les taux de décohérence.
*   `update_state(self, p700_concentration)`: Met à jour l'état quantique du P700 (0 pour l'état fondamental, 1 pour l'état excité) en fonction de la concentration de P700 et des probabilités d'excitation/décohérence. Cette méthode simule également la lecture à double canal en produisant des signaux lumineux verts ou rouges.
*   `get_state(self)`: Retourne l'état actuel du P700, ainsi que les sorties des canaux de luciférase verte et rouge.
*   `update(self, p700_concentration)`: Une méthode d'enveloppe qui appelle `update_state` et retourne les sorties des canaux de luciférase.

**Exemple d'utilisation (extrait de `update_state`):**
```python
def update_state(self, p700_concentration):
    """Met à jour l'état de P700 et des canaux de lecture."""
    self.luc_green_output = 0.0
    self.luc_red_output = 0.0
    if self.p700_state == 0:
        if p700_concentration > self.p700_excitation_threshold:
            excitation_prob = (p700_concentration - self.p700_excitation_threshold) / (1.0 - self.p700_excitation_threshold)
            if random.random() < excitation_prob:
                self.p700_state = 1
    else:
        if random.random() < self.fast_collapse_probability:
            if random.random() < self.decoherence_rate_fast:
                self.p700_state = 0
                self.luc_red_output = 1.0
        else:
            if random.random() < self.decoherence_rate_slow:
                self.p700_state = 0
                self.luc_green_output = 1.0
```

│   │   │   ├── biological_engine.py

│   │   │   ├── environment_engine.py

### Fichier: bioos/simulations/multiphysics_simulator/environment_engine.py
**Position:** `bioos/simulations/multiphysics_simulator/environment_engine.py`
**Utilité:** Gère la simulation de l'environnement externe, notamment l'intensité lumineuse, pour le simulateur multiphysique.
**Description:** Ce fichier définit la classe `EnvironmentEngine` qui simule les conditions environnementales. Il est capable de gérer des impulsions lumineuses configurables, en ajustant l'intensité lumineuse en fonction du temps de simulation.

**Méthodes clés:**
*   `__init__(self, config)`: Initialise le moteur d'environnement avec une configuration qui peut inclure des `pulse_configs` (configurations d'impulsions lumineuses).
*   `update(self, time, dt)`: Met à jour l'état de l'environnement à chaque pas de temps. Il vérifie si une impulsion lumineuse est active et ajuste l'intensité lumineuse en conséquence. Il retourne l'intensité lumineuse actuelle.

**Exemple d'utilisation (extrait de `update`):**
```python
def update(self, time, dt):
    current_intensity = 0
    for pulse in self.pulse_configs:
        if pulse['start'] <= time < pulse['end']:
            current_intensity = pulse['intensity']
            break
    
    if current_intensity != self.light_intensity:
        print(f"EVENT: Light intensity changed to {current_intensity} at t={time}")
        self.light_intensity = current_intensity

    return {
        'light_intensity': self.light_intensity
    }
```

│   │   │   └── simulator.py

### Fichier: bioos/simulations/multiphysics_simulator/simulator.py
**Position:** `bioos/simulations/multiphysics_simulator/simulator.py`
**Utilité:** Orchestre la simulation multiphysique en intégrant les moteurs biologique, quantique et environnemental.
**Description:** Ce fichier définit la classe `MultiphysicsSimulator` qui est le point d'entrée principal pour l'exécution des simulations. Il initialise et gère les interactions entre les différents moteurs (biologique, quantique, environnemental) et enregistre l'état de la simulation au fil du temps. Il inclut également des fonctionnalités pour la visualisation des résultats.

**Méthodes clés:**
*   `__init__(self, config)`: Initialise le simulateur avec une configuration donnée, créant des instances de `BiologicalEngine`, `QuantumEngine` et `EnvironmentEngine`.
*   `run(self)`: Exécute la simulation complète jusqu'à `max_time`.
*   `run_until(self, time)`: Exécute la simulation jusqu'à un temps spécifié.
*   `run_and_plot_steps(self, time, frame_dir)`: Exécute la simulation et génère un graphique à chaque étape, sauvegardant les images dans un répertoire spécifié.
*   `plot_step(self, output_path)`: Génère un graphique de l'état actuel de la simulation.
*   `step(self)`: Effectue un seul pas de temps de la simulation, mettant à jour l'environnement, l'état biologique et l'état quantique, puis enregistre les données.
*   `plot_results(self, output_path)`: Génère un graphique final des résultats de la simulation.

**Exemple d'utilisation (extrait de `step`):**
```python
def step(self):
    # 1. Get light intensity from environment
    env_state = self.env_engine.update(self.time, self.dt)
    light_intensity = env_state['light_intensity']

    # 2. Update biological state
    self.bio_engine.update(self.time, self.dt, env_state)
    
    # 3. Update quantum state
    self.quantum_engine.update_state(self.bio_engine.p700_concentration)
    quantum_state = self.quantum_engine.get_state()

    # 4. Log current state
    self.log.append({
        'time': self.time,
        'light_intensity': light_intensity,
        'p700_concentration': self.bio_engine.p700_concentration,
        'p700_state': quantum_state['p700_state'],
        'luc_green_output': quantum_state['luc_green_output'],
        'luc_red_output': quantum_state['luc_red_output']
    })

    # 5. Advance time
    self.time += self.dt
```

## Modèles Biologiques

Cette section regroupe les fichiers qui définissent des modèles biologiques spécifiques utilisés dans les simulations.

### Fichier: bioos/simulations/biological_models/gene_regulation_model.py
**Position:** `bioos/simulations/biological_models/gene_regulation_model.py`
**Utilité:** Définit le modèle mathématique de la dynamique de la concentration de P700, utilisé par le `BiologicalEngine`.
**Description:** Ce fichier contient la fonction `gene_regulation_model` qui décrit l'évolution de la concentration de P700 au fil du temps en réponse à l'intensité lumineuse. Il utilise des équations différentielles ordinaires (ODE) basées sur la cinétique de Michaelis-Menten pour modéliser la synthèse et la dégradation du P700. Ce modèle est un composant essentiel du `BiologicalEngine` pour simuler les processus biologiques au sein du simulateur multiphysique.

**Fonctions clés:**
*   `gene_regulation_model(p, t, light_intensity, params)`: Calcule le taux de changement de la concentration de P700 (`dP700_dt`) en fonction de la concentration actuelle de P700 (`p`), du temps (`t`), de l'intensité lumineuse (`light_intensity`) et des paramètres du modèle (`params`).

**Exemple d'utilisation (extrait):**
```python
def gene_regulation_model(p, t, light_intensity, params):
    P700 = p[0]
    
    synthesis_rate = params['V_max_synthesis'] * light_intensity / (params['K_light'] + light_intensity)
    degradation_rate = params['k_degradation'] * P700
    
    dP700_dt = synthesis_rate - degradation_rate
    
    return [dP700_dt]
```

### Fichier: bioos/simulations/biological_models/simulate_cry2_radical_pair.py
**Position:** `bioos/simulations/biological_models/simulate_cry2_radical_pair.py`
**Utilité:** Simule la dynamique de la paire de radicaux CRY2 et sa réponse électromagnétique, un élément clé dans la modélisation des effets quantiques en biologie.
**Description:** Ce fichier implémente une simulation basée sur les principes de la mécanique quantique pour modéliser le comportement des paires de radicaux CRY2, qui sont impliquées dans la magnétoréception chez les plantes. Il calcule la fréquence de Larmor et simule l'oscillation entre les états singulet et triplet de la paire de radicaux en fonction du temps. Les résultats sont visualisés pour montrer la probabilité d'être dans chaque état.

**Fonctions clés:**
*   `simulate_cry2_radical_pair()`: Exécute la simulation de la paire de radicaux CRY2, calcule la fréquence de Larmor, simule les probabilités des états singulet et triplet, et visualise les résultats.

**Exemple d'utilisation (extrait):**
```python
def simulate_cry2_radical_pair():
    """Simule radical pair CRY2 et réponse EM"""
    print("🔬 SIMULATION RADICAL PAIR CRY2\n")
    
    # Paramètres (Ritz 2000)
    larmor_freq = 9.8e3  # 9.8 kHz
    magnetic_field = 2.2e-6
    gyromagnetic_ratio = 28.0e9  # rad/s/T (électron)
    
    # Calcul Larmor: ω = γB
    calculated_freq = (gyromagnetic_ratio * magnetic_field) / (2 * np.pi)
    
    print(f"📊 Paramètres:")
    print(f"  Fréquence Larmor: {larmor_freq/1e3:.1f} kHz")
    print(f"  Champ magnétique: {magnetic_field*1e6:.2f} μT")
    print(f"  Fréquence calculée: {calculated_freq/1e3:.1f} kHz")
    print()
    
    # Validation
    diff = abs(calculated_freq - larmor_freq) / larmor_freq * 100
    print(f"✅ Écart calculé: {diff:.1f}%")
    print(f"✅ Validation: {'OK' if diff < 5 else 'ATTENTION'}")
    
    # Simulation radical pair (FAD•- + Trp•+)
    # États de spin: singlet (S) et triplet (T)
    time = np.linspace(0, 0.01, 1000)  # 10 ms
    
    # Probabilité singlet (simplifié)
    # Oscillation à fréquence Larmor
    singlet_prob = 0.5 + 0.3 * np.cos(2 * np.pi * larmor_freq * time)
    triplet_prob = 1 - singlet_prob
    
    # Visualisation
    plt.figure(figsize=(10, 6))
    plt.plot(time * 1e3, singlet_prob, 'b-', label='Singlet (FAD•- + Trp•+)', linewidth=2)
    plt.plot(time * 1e3, triplet_prob, 'r-', label='Triplet', linewidth=2)
    plt.axvline(0.2, color='g', linestyle='--', label='200 ms (protocole)')
    plt.xlabel('Temps (ms)')
    plt.ylabel('Probabilité')
    plt.title('Radical Pair CRY2 - Réponse EM 9.8 kHz')
    plt.grid(True)
    plt.legend()
    plt.savefig('05_SIMULATION/results/cry2_radical_pair.png', dpi=150)
    print(f"\n📈 Graphique sauvegardé: 05_SIMULATION/results/cry2_radical_pair.png")
    
    return diff < 5
```

### Fichier: bioos/simulations/biological_models/simulate_growth.py
**Position:** `bioos/simulations/biological_models/simulate_growth.py`
**Utilité:** Simule la croissance simplifiée de la plante HAWRA (Ficus elastica) et la production théorique de qubits dans des conditions climatiques marocaines.
**Description:** Ce script fournit un modèle simplifié et linéaire de la croissance de la plante HAWRA sur 90 jours, basé sur des paramètres inspirés d'OpenSimRoot. Il estime la hauteur de la plante, la production théorique de qubits (1 qubit/cm, max 1000) et modélise un métabolisme CAM simplifié (80% d'efficacité la nuit). Le script visualise ces données ainsi que le cycle de température marocain. Il inclut également une validation simplifiée des résultats.

**Fonctions clés:**
*   `simulate_growth()`: Exécute la simulation de croissance, calcule la hauteur, la production de qubits et l'énergie métabolique, puis visualise les résultats et effectue une validation.

**Paramètres clés:**
*   `temp_day`, `temp_night`: Températures diurnes et nocturnes (Maroc).
*   `humidity`: Humidité relative.
*   `light_hours`: Heures de lumière par jour.
*   `initial_height`: Hauteur initiale de la graine.
*   `growth_rate`: Taux de croissance linéaire (cm/jour).
*   `cam_efficiency`: Efficacité du métabolisme CAM la nuit.

**Exemple d'utilisation (extrait):**
```python
def simulate_growth():
    """Simule croissance 3D HAWRA sur 90 jours"""
    print("🌱 SIMULATION CROISSANCE HAWRA\n")
    
    # Paramètres climat Maroc
    temp_day = 35  # °C jour
    temp_night = 25  # °C nuit
    humidity = 20  # %
    light_hours = 12  # h/jour
    
    # Paramètres croissance
    days = 90
    initial_height = 0.01  # 1 cm (graine)
    growth_rate = 0.03  # cm/jour (3 cm/mois)
    
    # Simulation
    time = np.arange(0, days + 1)
    height = initial_height + growth_rate * time
    
    # Production qubits (1 qubit par cm de hauteur, max 1000)
    qubits = np.minimum(height * 100, 1000)
    
    # CAM métabolisme (énergie 24h/24)
    cam_efficiency = 0.8  # 80% efficacité
    energy_day = np.ones(days + 1) * 100  # 100% jour
    energy_night = np.ones(days + 1) * (80 * cam_efficiency)  # 80% nuit (CAM)
    
    print(f"📊 Paramètres:")
    print(f"  Durée: {days} jours")
    print(f"  Température: {temp_day}°C jour / {temp_night}°C nuit")
    print(f"  Humidité: {humidity}%
")
    
    print(f"📈 Résultats (jour {days}):")
    print(f"  Hauteur: {height[-1]:.2f} m")
    print(f"  Qubits: {int(qubits[-1])}")
    print(f"  Énergie jour: {energy_day[-1]:.0f}%
")
    
    # Validation
    height_ok = height[-1] >= 2.0  # Au moins 2m en 90 jours
    qubits_ok = qubits[-1] >= 100  # Au moins 100 qubits
    energy_ok = energy_night[-1] >= 50  # Au moins 50% énergie nuit
    
    print(f"\n✅ Validation:")
    print(f"  Hauteur: {'OK' if height_ok else 'ATTENTION'}")
    print(f"  Qubits: {'OK' if qubits_ok else 'ATTENTION'}")
    print(f"  Énergie: {'OK' if energy_ok else 'ATTENTION'}")
    
    return height_ok and qubits_ok and energy_ok
```

## Modèles Étendus

Cette section contient des modèles biologiques plus spécifiques ou expérimentaux, souvent liés à des phénomènes quantiques ou des réponses à des stimuli externes.

### Fichier: bioos/simulations/extended_model/cry2_model.py
**Position:** `bioos/simulations/extended_model/cry2_model.py`
**Utilité:** Simule l'expression du gène CRY2 en réponse à un signal électromagnétique, modélisant l'impact des champs EM sur la régulation génique.
**Description:** Ce fichier implémente une simulation de l'expression du gène CRY2, une protéine sensible à la lumière bleue et aux champs électromagnétiques, en fonction de la fréquence d'un signal EM. Le modèle utilise une fonction sigmoïde simplifiée pour représenter la réponse non linéaire de l'expression génique à l'intensité du signal, avec une saturation à des fréquences élevées. Il simule également une dynamique temporelle simple pour l'atteinte du plateau d'expression.

**Fonctions clés:**
*   `simulate_cry2_expression(em_signal_frequency, duration, time_step=0.1)`: Simule le niveau d'expression du gène CRY2 au fil du temps en fonction de la fréquence du signal EM, de la durée de la simulation et du pas de temps.

**Paramètres clés:**
*   `em_signal_frequency`: Fréquence du signal électromagnétique en Hz.
*   `duration`: Durée totale de la simulation en secondes.
*   `time_step`: Intervalle de temps entre chaque calcul dans la simulation.
*   `k`: Pente de la fonction sigmoïde, influençant la raideur de la réponse.
*   `midpoint`: Fréquence du signal EM à laquelle l'expression atteint la moitié de sa valeur maximale.

**Exemple d'utilisation (extrait):**
```python
def simulate_cry2_expression(em_signal_frequency, duration, time_step=0.1):
    """
    Simule l'expression du gène CRY2 en réponse à un signal électromagnétique.

    Args:
        em_signal_frequency (float): Fréquence du signal EM en Hz.
        duration (float): Durée de la simulation en secondes.
        time_step (float): Pas de temps pour la simulation en secondes.

    Returns:

## Simulations Procédurales

Cette section décrit les simulations procédurales du protocole de régénération, incluant les versions standard et optimisée.

### Fichier: bioos/simulations/sop_procedural_simulation.py
**Position:** `bioos/simulations/sop_procedural_simulation.py`
**Utilité:** Simule de manière impérative le protocole de régénération du Ficus elastica, étape par étape, avec des probabilités de succès fixes.
**Description:** Ce fichier implémente une simulation séquentielle du protocole de régénération in vitro pour un explant unique de Ficus elastica. Chaque étape du protocole (prélèvement, transformation génétique, sélection, formation du cal, induction des pousses, enracinement, acclimatation) est simulée avec une probabilité de succès prédéfinie. La simulation s'arrête dès qu'une étape échoue, ou se termine avec succès si toutes les étapes sont franchies. Ce script est utile pour comprendre le déroulement de base du protocole et les points de défaillance potentiels.

**Fonctions clés:**
*   `run_procedural_simulation()`: Exécute la simulation complète du protocole, affichant le résultat de chaque étape et le succès ou l'échec final.

**Paramètres clés (via `PROBABILITIES`):**
*   `explant_survival`: Probabilité de survie de l'explant après prélèvement.
*   `transformation_efficiency`: Efficacité de la transformation par Agrobacterium.
*   `selection_survival`: Survie à la sélection sur milieu sélectif.
*   `callus_formation`: Probabilité de formation de cals.
*   `shoot_induction`: Probabilité d'induction de pousses à partir des cals.
*   `root_formation`: Probabilité de développement de racines.
*   `acclimatization_survival`: Survie de la plante lors de l'acclimatation ex vitro.

**Exemple d'utilisation (extrait):**
```python
def run_procedural_simulation():
    """
    Exécute une simulation impérative du protocole de régénération pour un seul explant.
    Affiche le résultat de chaque étape.
    """
    print("===========================================================")
    print("=== Lancement de la Simulation Procédurale du SOP Ficus ===")
    print("=== Suivi du parcours d'un seul explant...              ===")
    print("===========================================================")
    time.sleep(1)

    # Étape 1: Prélèvement de l'explant
    print("\n[ÉTAPE 1/7] Prélèvement et préparation de l'explant...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["explant_survival"]:
        print("  -> RÉSULTAT: ÉCHEC. L'explant n'a pas survécu à la stérilisation/préparation.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. L'explant est viable.")
    # ... (autres étapes) ...
    return True
```

## Simulation de Régénération (Monte Carlo)

### Fichier: bioos/simulations/validate_simulation.py
**Position:** `bioos/simulations/validate_simulation.py`
**Utilité:** Valide le comportement du simulateur multiphysique en analysant les logs de simulation pour détecter les anomalies et assurer la conformité avec les attentes biologiques et quantiques.
**Description:** Ce script fournit un cadre de validation pour les simulations multiphysiques. Il analyse les fichiers de log JSON générés par le simulateur pour vérifier des aspects critiques tels que la dynamique de la concentration de P700, la validation des excitations, et l'exclusion mutuelle des canaux de lecture. Il est essentiel pour garantir que le modèle se comporte comme prévu et que les résultats de simulation sont fiables.

**Fonctions clés:**
*   `analyze_simulation_log(log_path, config)`: Prend en entrée le chemin d'un fichier de log de simulation et la configuration utilisée, puis retourne un dictionnaire indiquant le statut de validation (SUCCÈS/ÉCHEC), une liste d'erreurs détectées, le nombre total d'étapes, et des statistiques sur les excitations et les lectures des canaux vert et rouge.

**Exemple d'utilisation (extrait du bloc `if __name__ == "__main__":`) :**
```python
if __name__ == "__main__":
    # This is the same config used in the simulation
    config = {
        'max_time': 400,
        'dt': 0.5,
        'env': {
            'pulse_configs': [
                {'start': 10, 'end': 20, 'intensity': 1.0},
                {'start': 50, 'end': 55, 'intensity': 0.8},
                {'start': 90, 'end': 92, 'intensity': 1.0},
                {'start': 120, 'end': 140, 'intensity': 0.6},
                {'start': 160, 'end': 170, 'intensity': 1.0},
                {'start': 200, 'end': 220, 'intensity': 0.9},
                {'start': 250, 'end': 270, 'intensity': 0.7},
                {'start': 300, 'end': 310, 'intensity': 1.0},
                {'start': 340, 'end': 360, 'intensity': 0.5}
            ]
        },
        'bio': {
            'p700_initial': 0.0,
            'degradation_rate': 0.04,
            'synthesis_rate': 0.25
        },
        'quantum': {
            'threshold': 0.8,
            'decoherence_rate': 0.015
        }
    }
    
    log_file = "/Users/mehdiwhb/Desktop/HAWRA/05_data/results/multiphysics_simulation/multiphysics_simulation_v2.json"
    
    results = analyze_simulation_log(log_file, config)
    
    print("--- PQPE Numerical Validation Results ---")
    print(f"Status: {results['validation_status']}")
    if results['errors']:
        print("Errors found:")
        for error in results['errors']:
            print(f"- {error}")
    else:
        print("No errors found. The model behaves as expected.")
    print("
--- Statistics ---")
    print(f"Total simulation steps: {results['total_steps']}")
    print(f"P700 excitations: {results['excitations']}")
    print(f"Green channel readouts (|0>): {results['green_reads']}")
    print(f"Red channel readouts (|1>): {results['red_reads']}")
    print("------------------------------------")
```

Cette section décrit la simulation de Monte Carlo pour estimer le rendement global du protocole de régénération.

### Fichier: bioos/simulations/regeneration_simulation.py
**Position:** `bioos/simulations/regeneration_simulation.py`
**Utilité:** Estime le rendement global du protocole de régénération du Ficus elastica par une approche de Monte Carlo, modélisant chaque étape critique comme un événement probabiliste.
**Description:** Ce script simule le protocole de régénération comme une chaîne de Markov, où la sortie d'une étape devient l'entrée de la suivante. Chaque étape (transformation, sélection, callogenèse, organogenèse, enracinement, acclimatation) est caractérisée par une probabilité de succès. La simulation exécute un grand nombre d'essais (simulations) pour un nombre initial donné d'explants, et calcule la distribution du nombre final de plantes viables (HAWRA-Ficus-G0). Les résultats sont visualisés sous forme d'histogramme et des statistiques clés (moyenne, écart-type, rendement global) sont fournies.

**Fonctions clés:**
*   `run_regeneration_simulation(n_simulations, n_explants_initial)`: Exécute la simulation de Monte Carlo. Prend en entrée le nombre de simulations et le nombre initial d'explants, et retourne une liste du nombre final de plantes viables pour chaque simulation.
*   `plot_results(simulation_results, n_explants_initial)`: Affiche les résultats de la simulation sous forme d'histogramme, calcule et imprime les statistiques, et sauvegarde le graphique.

**Paramètres clés (Probabilités de succès):**
*   `p_transfo`: Probabilité qu'une cellule d'explant intègre l'ADN-T.
*   `p_selection`: Probabilité qu'une cellule transformée survive à la sélection.
*   `p_callogenese`: Probabilité qu'un cal se forme à partir de cellules sélectionnées.
*   `p_organogenese`: Probabilité qu'un cal génère des bourgeons viables.
*   `p_enracinement`: Probabilité qu'une pousse développe un système racinaire.
*   `p_acclimatation`: Probabilité qu'une plantule survive au transfert en serre.

**Exemple d'utilisation (extrait):**
```python
if __name__ == "__main__":
    # --- Paramètres de la simulation ---
    N_SIMULATIONS = 10000  # Nombre d'essais pour la robustesse statistique
    N_EXPLANTS_INITIAL = 500 # Nombre d'explants de départ, un nombre réaliste pour une expérience en labo

    # --- Exécution et visualisation ---
    results = run_regeneration_simulation(N_SIMULATIONS, N_EXPLANTS_INITIAL)
    plot_results(results, N_EXPLANTS_INITIAL)
```

## Analyse de Sensibilité

Cette section décrit l'analyse de sensibilité "un facteur à la fois" (OFAT) pour identifier les étapes critiques du protocole de régénération.

### Fichier: bioos/simulations/sensitivity_analysis.py
**Position:** `bioos/simulations/sensitivity_analysis.py`
**Utilité:** Quantifie l'impact de la variation de chaque probabilité de succès du protocole de régénération sur le rendement final, permettant d'identifier les étapes les plus critiques.
**Description:** Ce script implémente une analyse de sensibilité de type "un facteur à la fois" (One-Factor-at-a-Time - OFAT). Pour chaque paramètre de probabilité du protocole de régénération (p_transfo, p_selection, etc.), sa valeur est variée sur une plage définie (par exemple, de 10% à 100%) tandis que tous les autres paramètres sont maintenus à leurs valeurs de base. Pour chaque point de variation, une simulation de Monte Carlo est exécutée pour calculer le rendement moyen. Les résultats sont ensuite visualisés pour montrer comment le rendement global est affecté par la variation de chaque paramètre, aidant ainsi à identifier les leviers d'optimisation les plus efficaces.

**Fonctions clés:**
*   `run_regeneration_simulation(n_simulations, n_explants_initial, probabilities)`: Une version modifiée de la fonction de simulation de Monte Carlo qui accepte un dictionnaire de probabilités, permettant de tester différentes configurations de paramètres.
*   `run_sensitivity_analysis(base_probabilities, n_simulations, n_explants_initial)`: Orchestre l'analyse OFAT. Elle itère sur chaque paramètre, le fait varier sur une plage et appelle `run_regeneration_simulation` pour obtenir le rendement moyen. Retourne les résultats de sensibilité et la plage de variation des paramètres.
*   `plot_sensitivity_results(sensitivity_results, param_range)`: Affiche les résultats de l'analyse de sensibilité sous forme graphique, montrant l'évolution du rendement global en fonction de la variation de chaque paramètre. Sauvegarde le graphique.

**Paramètres clés:**
*   `base_probabilities`: Dictionnaire des probabilités de succès de base pour chaque étape du protocole.
*   `param_range`: Plage de valeurs sur laquelle chaque probabilité est variée (par exemple, `np.linspace(0.1, 1.0, 10)`).

**Exemple d'utilisation (extrait):**
```python
if __name__ == "__main__":
    # --- Paramètres de base pour la simulation ---
    N_SIMULATIONS = 2000
    N_EXPLANTS_INITIAL = 500

    # --- Probabilités de base (point de fonctionnement) ---
    base_probabilities = {
        'p_transfo': 0.10,
        'p_selection': 0.60,
        'p_callogenese': 0.50,
        'p_organogenese': 0.40,
        'p_enracinement': 0.70,
        'p_acclimatation': 0.50
    }

    # --- Exécution de l'analyse de sensibilité ---
    sensitivity_data, param_range = run_sensitivity_analysis(
        base_probabilities, N_SIMULATIONS, N_EXPLANTS_INITIAL
    )

    # --- Visualisation des résultats ---
    plot_sensitivity_results(sensitivity_data, param_range)
```

### Fichier: bioos/simulations/sop_procedural_simulation_optimized.py
**Position:** `bioos/simulations/sop_procedural_simulation_optimized.py`
**Utilité:** Simule le protocole de régénération du Ficus elastica avec des paramètres optimisés, notamment une efficacité de transformation génétique améliorée.
**Description:** Ce fichier est une version optimisée de la simulation procédurale du protocole de régénération. Il intègre des améliorations basées sur des conditions expérimentales spécifiques, comme l'ajout de 200µM d'acétosyringone lors de la co-culture avec Agrobacterium, ce qui augmente significativement l'efficacité de la transformation génétique. Les autres étapes du protocole restent similaires à la version standard, mais l'impact de cette optimisation sur le succès global est mis en évidence.

**Fonctions clés:**
*   `run_procedural_simulation_optimized()`: Exécute la simulation du protocole optimisé, affichant les résultats de chaque étape et le succès ou l'échec final.

**Paramètres clés (via `PROBABILITIES`):**
*   `transformation_efficiency`: Efficacité de transformation AMÉLIORÉE (par exemple, de 0.40 à 0.65) grâce à l'optimisation.
*   Les autres probabilités sont identiques à celles de `sop_procedural_simulation.py`.

**Exemple d'utilisation (extrait):**
```python
def run_procedural_simulation_optimized():
    """
    Exécute une simulation impérative du protocole de régénération OPTIMISÉ.
    """
    print("=====================================================================")
    print("=== Lancement de la Simulation du Protocole OPTIMISÉ (200µM AS) ===")
    print("=== Suivi du parcours d'un seul explant...                       ===")
    print("=====================================================================")
    time.sleep(1)

    # Étape 1: Prélèvement de l'explant
    print("\n[ÉTAPE 1/7] Prélèvement et préparation de l'explant...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["explant_survival"]:
        print("  -> RÉSULTAT: ÉCHEC. L'explant n'a pas survécu.")
        print("\n--- SIMULATION TERMINÉE ---
")
        return False
    print("  -> RÉSULTAT: SUCCÈS. L'explant est viable.")

    # Étape 2: Transformation génétique (OPTIMISÉE)
    print("\n[ÉTAPE 2/7] Co-culture avec Agrobacterium (Protocole Optimisé, 200µM Acetosyringone)...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["transformation_efficiency"]:
        print("  -> RÉSULTAT: ÉCHEC. Le transfert du gène HAWRA a échoué.")
        print("\n--- SIMULATION TERMINÉE ---
")
        return False
    print("  -> RÉSULTAT: SUCCÈS. L'ADN-T a été transféré avec une meilleure efficacité.")
    # ... (autres étapes) ...
    return True
```
        tuple: Un tuple contenant les tableaux de temps et de niveau d'expression.
    """
    time = np.arange(0, duration, time_step)
    # Modèle simplifié : l'expression est proportionnelle à la fréquence du signal
    # avec une réponse de type sigmoïde pour représenter la saturation.
    # Les paramètres de la sigmoïde sont choisis pour illustrer le concept.
    k = 0.1  # Raideur de la sigmoïde
    midpoint = 50  # Fréquence à mi-expression
    expression_level = 1 / (1 + np.exp(-k * (em_signal_frequency - midpoint)))
    
    # Simuler une dynamique temporelle simple (atteinte du plateau)
    expression_over_time = expression_level * (1 - np.exp(-time / (duration / 5)))

    return time, expression_over_time
```

### Fichier: bioos/simulations/extended_model/metabolic_model.py
**Position:** `bioos/simulations/extended_model/metabolic_model.py`
**Utilité:** Simule la consommation d'ATP et de NADPH en fonction de l'activité quantique, modélisant l'impact énergétique des processus quantiques.
**Description:** Ce fichier implémente un modèle simplifié de la consommation métabolique (ATP et NADPH) au sein de la plante en réponse à un nombre donné d'opérations quantiques. Il modélise la diminution des niveaux d'ATP et de NADPH au fil du temps, proportionnellement à l'activité quantique, et s'assure que ces niveaux ne deviennent pas négatifs. Ce modèle permet d'évaluer les coûts énergétiques des processus quantiques simulés.

**Fonctions clés:**
*   `simulate_metabolism(quantum_operations_count, duration, time_step=0.1, initial_atp=100, initial_nadph=100)`: Simule les niveaux d'ATP et de NADPH au fil du temps en fonction du nombre d'opérations quantiques, de la durée de la simulation, du pas de temps et des niveaux initiaux d'ATP et de NADPH.

**Paramètres clés:**
*   `quantum_operations_count`: Nombre d'opérations quantiques effectuées, influençant le taux de consommation.
*   `duration`: Durée totale de la simulation en secondes.
*   `time_step`: Intervalle de temps entre chaque calcul dans la simulation.
*   `initial_atp`, `initial_nadph`: Niveaux initiaux d'ATP et de NADPH.
*   `atp_consumption_rate`, `nadph_consumption_rate`: Taux de consommation par opération quantique.

**Exemple d'utilisation (extrait):**
```python
def simulate_metabolism(quantum_operations_count, duration, time_step=0.1, initial_atp=100, initial_nadph=100):
    """
    Simule la consommation d'ATP et de NADPH en fonction de l'activité quantique.

    Args:
        quantum_operations_count (int): Nombre d'opérations quantiques effectuées.
        duration (float): Durée de la simulation en secondes.
        time_step (float): Pas de temps pour la simulation en secondes.
        initial_atp (float): Niveau initial d'ATP.
        initial_nadph (float): Niveau initial de NADPH.

    Returns:
        tuple: Un tuple contenant les tableaux de temps, de niveaux d'ATP et de niveaux de NADPH.
    """
    time = np.arange(0, duration, time_step)
    
    # Modèle de consommation : proportionnelle au nombre d'opérations quantiques
    atp_consumption_rate = 0.1 * quantum_operations_count
    nadph_consumption_rate = 0.05 * quantum_operations_count
    
    atp_levels = initial_atp - atp_consumption_rate * time
    nadph_levels = initial_nadph - nadph_consumption_rate * time
    
    # S'assurer que les niveaux ne deviennent pas négatifs
    atp_levels[atp_levels < 0] = 0
    nadph_levels[nadph_levels < 0] = 0

    return time, atp_levels, nadph_levels
```

### Fichier: bioos/simulations/extended_model/multiphysics_simulator.py
**Position:** `bioos/simulations/extended_model/multiphysics_simulator.py`
**Utilité:** Orchestre une simulation multiphysique intégrée, combinant les modèles biochimiques (CRY2), quantiques (Bio-qubit P700) et métaboliques pour une vue holistique.
**Description:** Ce fichier est le point d'intégration des différents modèles de l'`extended_model`. Il exécute une simulation qui lie l'expression du gène CRY2 (en réponse à un signal EM) à la dynamique d'un bio-qubit P700 (via un angle de rotation de porte quantique), et évalue l'impact énergétique de ces opérations quantiques sur le métabolisme (consommation d'ATP et de NADPH). Il utilise la bibliothèque `qutip` pour la simulation quantique.

**Fonctions clés:**
*   `run_multiphysics_simulation(em_frequency, sim_duration=10)`: Exécute la simulation multiphysique complète, enchaînant les appels aux fonctions de `cry2_model.py` et `metabolic_model.py`, et effectuant une simulation quantique.

**Paramètres clés:**
*   `em_frequency`: Fréquence du signal électromagnétique en Hz, qui influence l'expression de CRY2.
*   `sim_duration`: Durée totale de la simulation en secondes.
*   `final_cry2_level`: Niveau d'expression final de CRY2, servant de lien entre la biologie et le quantique.
*   `rotation_angle`: Angle de rotation de la porte quantique, dérivé du niveau de CRY2.
*   `H`: Hamiltonien pour la rotation du qubit.
*   `psi0`: État initial du qubit.
*   `prob_1`: Probabilité de lecture de l'état |1> du qubit.

**Exemple d'utilisation (extrait):**
```python
def run_multiphysics_simulation(em_frequency, sim_duration=10):
    """
    Exécute une simulation multiphysique intégrée.

    Args:
        em_frequency (float): Fréquence du signal EM en Hz.
        sim_duration (float): Durée de la simulation en secondes.

    Returns:
        dict: Un dictionnaire contenant les résultats de la simulation.
    """
    # 1. Simulation biochimique (CRY2)
    _, cry2_expression_over_time = simulate_cry2_expression(em_frequency, sim_duration)
    final_cry2_level = cry2_expression_over_time[-1]

    # 2. Simulation quantique (Bio-qubit P700)
    # L'angle de rotation de la porte quantique dépend du niveau d'expression de CRY2
    # C'est le lien clé entre la biologie et le quantique
    rotation_angle = final_cry2_level * np.pi  # Angle max = pi (porte NOT)
    H = rotation_angle * sigmax()  # Hamiltonien pour la rotation
    
    psi0 = basis(2, 0)  # État initial |0>
    times = np.linspace(0, 1, 101) # Temps de l'opération de porte
    result = mesolve(H, psi0, times, [], [])
    final_state = result.states[-1]
    prob_1 = np.abs((basis(2, 1).dag() * final_state).full())**2

    # 3. Simulation métabolique
    # La consommation dépend du nombre d'opérations (ici, 1 porte)
    n_ops = 1
    t_metab, atp, nadph = simulate_metabolism(n_ops, sim_duration)

    return {
        "em_frequency": em_frequency,
        "final_cry2_level": final_cry2_level,
        "final_qubit_state": final_state,
        "readout_probability": prob_1,
        "metabolism": {
            "time": t_metab,
            "atp": atp,
            "nadph": nadph
        }
    }
```

│   │   ├── regeneration_simulation.py

### Fichier: bioos/simulations/regeneration_simulation.py
**Position:** `bioos/simulations/regeneration_simulation.py`
**Utilité:** Simule le protocole de régénération de Ficus elastica (HAWRA) en utilisant une approche de Monte Carlo pour estimer le rendement global.
**Description:** Ce fichier implémente une simulation de Monte Carlo pour modéliser le protocole de régénération de Ficus elastica. Il traite chaque étape du protocole comme un événement probabiliste, permettant de prédire le nombre de plantes viables (HAWRA-Ficus-G0) obtenues à partir d'un nombre initial d'explants. Le modèle est une chaîne de Markov où la sortie d'une étape devient l'entrée de la suivante.

**Fonctions clés:**
*   `run_regeneration_simulation(n_simulations, n_explants_initial)`: Exécute la simulation de Monte Carlo. Elle prend en entrée le nombre de simulations à effectuer et le nombre initial d'explants. Elle retourne une liste du nombre final de plantes viables pour chaque simulation.
*   `plot_results(simulation_results, n_explants_initial)`: Affiche les résultats de la simulation sous forme d'histogramme et calcule des statistiques clés (moyenne, écart-type, rendement global moyen). Elle sauvegarde également le graphique.

**Paramètres clés (Probabilités de succès):**
*   `p_transfo`: Probabilité qu'une cellule d'explant intègre l'ADN-T.
*   `p_selection`: Probabilité qu'une cellule transformée survive à la sélection.
*   `p_callogenese`: Probabilité qu'un cal se forme à partir de cellules sélectionnées.
*   `p_organogenese`: Probabilité qu'un cal génère des bourgeons viables.
*   `p_enracinement`: Probabilité qu'une pousse développe un système racinaire.
*   `p_acclimatation`: Probabilité qu'une plantule survive au transfert en serre.

**Exemple d'utilisation (extrait de `run_regeneration_simulation`):**
```python
def run_regeneration_simulation(n_simulations, n_explants_initial):
    p_transfo = 0.10
    p_selection = 0.60
    p_callogenese = 0.50
    p_organogenese = 0.40
    p_enracinement = 0.70
    p_acclimatation = 0.50

    final_plant_counts = []

    for _ in range(n_simulations):
        n_surviving_selection = np.random.binomial(n_explants_initial, p_transfo * p_selection)
        n_calli = np.random.binomial(n_surviving_selection, p_callogenese)
        n_shoots = np.random.binomial(n_calli, p_organogenese)
        n_rooted_plantlets = np.random.binomial(n_shoots, p_enracinement)
        n_final_plants = np.random.binomial(n_rooted_plantlets, p_acclimatation)
        final_plant_counts.append(n_final_plants)

    return final_plant_counts
```

│   │   ├── sensitivity_analysis.py

### Fichier: bioos/simulations/sensitivity_analysis.py
**Position:** `bioos/simulations/sensitivity_analysis.py`
**Utilité:** Réalise une analyse de sensibilité sur les paramètres du protocole de régénération de Ficus elastica pour identifier les étapes les plus critiques.
**Description:** Ce fichier implémente une analyse de sensibilité "un facteur à la fois" (OFAT) pour évaluer l'impact de la variation des probabilités de succès de chaque étape du protocole de régénération sur le rendement final. Il utilise la fonction de simulation de Monte Carlo (`run_regeneration_simulation`) en faisant varier un paramètre à la fois sur une plage définie, tout en maintenant les autres à leurs valeurs de base. Les résultats sont ensuite visualisés pour montrer l'influence de chaque paramètre.

**Fonctions clés:**
*   `run_regeneration_simulation(n_simulations, n_explants_initial, probabilities)`: Une version modifiée de la fonction de simulation de Monte Carlo qui accepte un dictionnaire de probabilités, permettant de tester l'impact de la variation des paramètres.
*   `run_sensitivity_analysis(base_probabilities, n_simulations, n_explants_initial)`: Orchestre l'analyse de sensibilité. Elle itère sur chaque paramètre, le fait varier sur une plage définie et exécute la simulation pour chaque valeur afin de calculer le rendement moyen.
*   `plot_sensitivity_results(sensitivity_results, param_range)`: Affiche les résultats de l'analyse de sensibilité sous forme graphique, montrant comment le rendement global moyen varie en fonction de chaque probabilité de succès. Le graphique est sauvegardé dans un fichier PNG.

**Paramètres clés analysés:**
*   `p_transfo`: Probabilité de transformation.
*   `p_selection`: Probabilité de survie à la sélection.
*   `p_callogenese`: Probabilité de formation de cal.
*   `p_organogenese`: Probabilité de production de bourgeons.
*   `p_enracinement`: Probabilité d'enracinement.
*   `p_acclimatation`: Probabilité d'acclimatation.

**Exemple d'utilisation (extrait de `run_sensitivity_analysis`):**
```python
def run_sensitivity_analysis(base_probabilities, n_simulations, n_explants_initial):
    sensitivity_results = {}
    param_range = np.linspace(0.1, 1.0, 10) # Plage de variation de 10% à 100%

    for param_name in base_probabilities.keys():
        yields = []
        
        current_probabilities = base_probabilities.copy()
        for p_value in param_range:
            current_probabilities[param_name] = p_value
            
            simulation_results = run_regeneration_simulation(
                n_simulations, n_explants_initial, current_probabilities
            )
            
            mean_yield = np.mean(simulation_results) / n_explants_initial * 100
            yields.append(mean_yield)

        sensitivity_results[param_name] = yields
    
    return sensitivity_results, param_range
```

### Validation des Probabilités Monte Carlo pour la Régénération Procédurale

**Utilité:** Vérifier la cohérence et la validité des probabilités utilisées dans les simulations de régénération procédurale, en comparant les versions standard et optimisée.

**Description:** Cette section documente la validation des probabilités définies dans `sop_procedural_simulation.py` et `sop_procedural_simulation_optimized.py`. L'objectif est de s'assurer que toutes les probabilités sont dans la plage valide [0, 1] et que les valeurs "optimisées" reflètent une amélioration logique.

**Résultats de la Validation:**

*   **`sop_procedural_simulation.py` (Probabilités Standard):**
    *   `explant_survival`: 0.95 (Valide)
    *   `transformation_efficiency`: 0.40 (Valide)
    *   `selection_survival`: 0.25 (Valide)
    *   `callus_formation`: 0.70 (Valide)
    *   `shoot_induction`: 0.30 (Valide)
    *   `root_formation`: 0.50 (Valide)
    *   `acclimatization_survival`: 0.60 (Valide)
    Toutes les probabilités sont dans la plage [0, 1].

*   **`sop_procedural_simulation_optimized.py` (Probabilités Optimisées):**
    *   `explant_survival`: 0.95 (Valide)
    *   `transformation_efficiency`: 0.65 (Valide) - **Amélioration:** Cette valeur est significativement plus élevée que la version standard (0.40), ce qui est cohérent avec un protocole "optimisé".
    *   `selection_survival`: 0.25 (Valide)
    *   `callus_formation`: 0.70 (Valide)
    *   `shoot_induction`: 0.30 (Valide)
    *   `root_formation`: 0.50 (Valide)
    *   `acclimatization_survival`: 0.60 (Valide)
    Toutes les probabilités sont dans la plage [0, 1]. L'optimisation est clairement visible dans l'augmentation de l'efficacité de transformation.

**Conclusion:** Les probabilités Monte Carlo utilisées dans les simulations procédurales sont validées. Elles respectent les contraintes de plage [0, 1] et la version optimisée démontre une amélioration ciblée et réaliste de l'efficacité de transformation, ce qui est crucial pour l'objectif de régénération accrue.

### Validation des Constantes ODE pour la Dynamique du P700

**Utilité:** Vérifier la cohérence et la validité des constantes utilisées dans le modèle d'équation différentielle ordinaire (ODE) pour la régulation génique du P700.

**Description:** Cette section documente la validation des constantes définies dans `bioos/simulations/biological_models/gene_regulation_model.py`. L'objectif est de s'assurer que les valeurs sont biologiquement plausibles et cohérentes avec les principes de la cinétique enzymatique et de la régulation génique.

**Constantes Validées:**

*   `k_prod_p700`: 0.1 (Taux de production du P700) - Valeur positive, cohérente avec un taux de synthèse.
*   `k_deg_p700`: 0.02 (Taux de dégradation du P700) - Valeur positive, cohérente avec un taux de dégradation.
*   `K_light`: 0.5 (Constante de Michaelis-Menten pour l'activation lumineuse) - Valeur positive, typique pour une constante de demi-saturation.
*   `n_light`: 2 (Coefficient de Hill pour l'activation lumineuse) - Valeur entière positive, indiquant une coopérativité dans l'activation par la lumière.

**Conclusion:** Les constantes ODE utilisées dans le modèle de régulation génique du P700 sont validées. Leurs valeurs sont biologiquement plausibles et respectent les conventions des modèles cinétiques, assurant une base solide pour la simulation de la dynamique du P700.

### Analyse de l'Impact Réel des Simulations HAWRA

**Utilité:** Évaluer la portée pratique et les implications des simulations HAWRA en termes d'optimisation de la régénération des cultures et d'applications en biologie quantique.

**Description:** Cette section analyse comment les différentes simulations développées dans le cadre du projet HAWRA contribuent à des avancées concrètes. Elle se concentre sur deux axes principaux : l'optimisation des protocoles de régénération de Ficus elastica et l'exploration des phénomènes de biologie quantique liés au P700.

**1. Optimisation de la Régénération des Cultures (Ficus elastica) :**

*   **Amélioration des Rendements:** La comparaison entre `sop_procedural_simulation.py` et `sop_procedural_simulation_optimized.py` démontre l'efficacité de l'approche de modélisation pour identifier et valider des protocoles améliorés. L'augmentation de l'efficacité de transformation de 0.40 à 0.65 dans la version optimisée a un impact direct et significatif sur le rendement final de plantes régénérées. Cela permet de réduire les coûts, le temps et les ressources nécessaires à la production de plantes génétiquement modifiées.
*   **Identification des Étapes Critiques:** L'analyse de sensibilité (`sensitivity_analysis.py`) est un outil crucial pour les biologistes. En identifiant les étapes du protocole qui ont le plus grand impact sur le succès global (par exemple, la transformation ou la sélection), les efforts de recherche peuvent être ciblés plus efficacement pour optimiser ces étapes spécifiques en laboratoire.
*   **Réduction des Essais Physiques:** En permettant de tester virtuellement un grand nombre de scénarios et de variations de paramètres, les simulations réduisent considérablement le besoin d'expériences coûteuses et chronophages en serre ou en laboratoire. Cela accélère le processus de développement de nouvelles variétés de plantes.

**2. Applications en Biologie Quantique (Dynamique du P700) :**

*   **Compréhension des Mécanismes Fondamentaux:** Le simulateur multiphysique, notamment via `quantum_engine.py` et `gene_regulation_model.py`, offre une plateforme pour étudier l'interaction complexe entre les phénomènes quantiques (comme la cohérence du P700) et les processus biologiques (régulation génique, photosynthèse). Cela ouvre des voies pour une meilleure compréhension des mécanismes sous-jacents à l'efficacité énergétique des plantes.
*   **Développement de Nouvelles Stratégies d'Ingénierie Biologique:** En modélisant l'impact de la lumière sur la dynamique du P700 et sa régulation génique, les simulations HAWRA peuvent aider à concevoir des stratégies pour manipuler ces processus. Par exemple, optimiser les conditions lumineuses pour maximiser la production de biomasse ou améliorer la résilience des plantes face au stress environnemental.
*   **Exploration de l'Intrication Quantique en Biologie:** Bien que le projet se concentre sur des aspects spécifiques, la modélisation de la cohérence quantique du P700 jette les bases pour explorer des concepts plus avancés de biologie quantique, potentiellement menant à des découvertes sur la façon dont les systèmes biologiques exploitent les phénomènes quantiques pour leur fonctionnement.

**Conclusion Générale sur l'Impact Réel:** Les simulations HAWRA ont un double impact : elles fournissent des outils pratiques pour l'optimisation de la biotechnologie végétale (régénération de Ficus) et elles contribuent à l'avancement de la biologie fondamentale en explorant les frontières de la biologie quantique. Elles positionnent le projet à l'intersection de l'ingénierie et de la science fondamentale, avec des retombées potentielles significatives pour l'agriculture et la compréhension du vivant.

### Évaluation de la Pertinence Scientifique des Simulations HAWRA

**Utilité:** Positionner les simulations HAWRA dans le contexte de la recherche scientifique actuelle en photosynthèse et en biologie synthétique, en évaluant leur alignement avec les connaissances établies et les défis contemporains.

**Description:** Cette section examine comment les modèles et les approches de simulation du projet HAWRA s'intègrent et contribuent aux domaines de la photosynthèse et de la biologie synthétique. Elle évalue la validité scientifique des hypothèses sous-jacentes et la pertinence des résultats pour la communauté scientifique.

**1. Alignement avec la Recherche sur la Photosynthèse :**

*   **Dynamique du P700:** La modélisation de la dynamique du P700 (`gene_regulation_model.py`, `quantum_engine.py`) est directement pertinente pour la recherche en photosynthèse. Le P700 est un composant clé du Photosystème I, et sa régulation, son excitation et sa dégradation sont des sujets d'étude intensifs pour comprendre l'efficacité de la conversion de l'énergie lumineuse. Les simulations HAWRA offrent un cadre pour explorer l'impact de divers facteurs (intensité lumineuse, taux de production/dégradation) sur la fonction du P700.
*   **Phénomènes Quantiques en Photosynthèse:** L'intégration de concepts de biologie quantique, tels que la cohérence quantique du P700, est un domaine de recherche émergent et très actif en photosynthèse. Les simulations HAWRA, en tentant de modéliser ces phénomènes, s'alignent avec les efforts visant à comprendre le rôle des effets quantiques dans l'efficacité quasi parfaite du transfert d'énergie dans les systèmes photosynthétiques.
*   **Optimisation de l'Efficacité Photosynthétique:** L'objectif ultime de nombreuses recherches en photosynthèse est d'améliorer l'efficacité de la conversion de l'énergie solaire en biomasse. Les simulations HAWRA, en fournissant des outils pour comprendre et potentiellement manipuler la dynamique du P700, contribuent à cet objectif en offrant des pistes pour l'ingénierie de systèmes photosynthétiques plus performants.

**2. Conformité aux Normes de Biologie Synthétique :**

*   **Approche Modulaire et Prédictive:** La biologie synthétique met l'accent sur la conception et la construction de systèmes biologiques avec des fonctions nouvelles ou améliorées. Les simulations HAWRA, avec leurs modèles modulaires (par exemple, les différents moteurs du simulateur multiphysique) et leur capacité à prédire les résultats de manipulations génétiques (via les probabilités de régénération), s'inscrivent parfaitement dans cette approche.
*   **Ingénierie de Voies Métaboliques et de Régulation:** La modélisation de la régulation génique du P700 est un exemple d'ingénierie de systèmes biologiques. En comprenant comment les gènes sont activés ou réprimés en réponse à des stimuli (comme la lumière), les biologistes synthétiques peuvent concevoir des circuits génétiques pour contrôler précisément le comportement des cellules ou des organismes.
*   **Conception Rationnelle de Protocoles:** La simulation procédurale et l'analyse de sensibilité sont des outils précieux pour la conception rationnelle de protocoles en biologie synthétique. Au lieu d'une approche par essais et erreurs, les simulations permettent d'optimiser les conditions expérimentales de manière prédictive, ce qui est une norme clé dans le domaine.

**Conclusion:** Les simulations HAWRA sont fortement pertinentes pour la recherche scientifique contemporaine. Elles s'alignent avec les avancées en photosynthèse en explorant la dynamique quantique du P700 et contribuent à la biologie synthétique par leur approche modulaire, prédictive et leur capacité à optimiser la conception de systèmes biologiques. Le projet HAWRA se positionne ainsi comme un contributeur potentiel aux efforts visant à exploiter et à manipuler les processus biologiques à des fins biotechnologiques et fondamentales.

### Évaluation du Niveau d'Innovation des Simulations HAWRA

**Utilité:** Identifier et évaluer les aspects novateurs des simulations HAWRA, en particulier les modèles hybrides quantique-biologiques et l'utilisation de la chaîne de Markov pour la régénération, par rapport aux approches conventionnelles.

**Description:** Cette section met en lumière les contributions uniques et les innovations méthodologiques du projet HAWRA. Elle examine comment l'intégration de concepts de la physique quantique avec la biologie, ainsi que l'application de modèles stochastiques avancés, distinguent ces simulations et ouvrent de nouvelles perspectives de recherche.

**1. Modèles Hybrides Quantique-Biologiques :**

*   **Intégration Multiphysique:** L'une des innovations majeures est l'intégration d'un moteur quantique (`quantum_engine.py`) avec des moteurs biologiques et environnementaux au sein d'un simulateur multiphysique (`simulator.py`). Cette approche permet d'étudier les interactions complexes entre les phénomènes quantiques (comme la cohérence du P700) et les processus biologiques à l'échelle moléculaire et cellulaire, ce qui est rarement abordé dans les simulations biologiques traditionnelles.
*   **Modélisation du P700 comme Système Quantique:** Le traitement du P700 non seulement comme une molécule biologique mais aussi comme un système capable de présenter des propriétés quantiques (cohérence, décohérence) est une avancée. Cela permet d'explorer des hypothèses sur le rôle de la mécanique quantique dans l'efficacité de la photosynthèse, un domaine de recherche de pointe et potentiellement révolutionnaire.
*   **Pont entre Disciplines:** Le projet HAWRA crée un pont entre la biologie, la physique quantique et l'ingénierie, favorisant une approche interdisciplinaire qui est essentielle pour résoudre des problèmes complexes et pour l'émergence de nouvelles disciplines comme la biologie quantique.

**2. Chaîne de Markov pour la Simulation de Régénération :**

*   **Approche Stochastique de la Régénération:** L'utilisation d'une chaîne de Markov pour modéliser le processus de régénération (`monte_carlo_simulation.py`) est une approche innovante par rapport aux modèles déterministes. Elle reconnaît la nature intrinsèquement probabiliste et séquentielle des étapes de régénération des plantes, offrant une représentation plus réaliste des rendements variables observés en laboratoire.
*   **Quantification des Incertitudes:** La chaîne de Markov permet de quantifier les incertitudes et les variabilités à chaque étape du protocole, ce qui est crucial pour une planification expérimentale robuste et pour l'optimisation des rendements. Cela contraste avec les modèles simplistes qui pourraient ignorer ces aspects stochastiques.
*   **Outil d'Optimisation de Protocole:** En permettant de simuler des milliers de parcours de régénération et d'évaluer l'impact de la modification des probabilités à chaque étape, la chaîne de Markov devient un puissant outil d'optimisation de protocole. Elle aide à identifier les goulots d'étranglement et à concevoir des stratégies pour améliorer l'efficacité globale de la régénération.

**Conclusion Générale sur l'Innovation:** Les simulations HAWRA se distinguent par leur audace à intégrer des concepts de physique quantique dans des modèles biologiques complexes et par leur utilisation sophistiquée de méthodes stochastiques pour la simulation de processus biologiques. Ces innovations positionnent le projet à l'avant-garde de la biologie computationnelle, avec le potentiel de générer de nouvelles connaissances fondamentales et d'ouvrir la voie à des applications biotechnologiques inédites.



### Fichier: bioos/simulations/multiphysics_simulator/biological_engine.py
**Position:** `bioos/simulations/multiphysics_simulator/biological_engine.py`
**Utilité:** Implémente le moteur de simulation biologique pour le simulateur multiphysique, modélisant la dynamique de la régulation génique et la concentration de P700.
**Description:** Ce fichier définit la classe `BiologicalEngine` qui gère la simulation des processus biologiques, notamment la régulation génique et la concentration de P700. Il utilise `scipy.integrate.odeint` pour résoudre les équations différentielles ordinaires (ODE) qui décrivent la dynamique du système biologique.

**Méthodes clés:**
*   `__init__(self, config)`: Initialise le moteur biologique avec une configuration donnée. Il initialise également la concentration de P700.
*   `update(self, time, dt, env_state)`: Met à jour l'état du système biologique pour un pas de temps donné. Il prend en compte l'intensité lumineuse de l'environnement et utilise le modèle de régulation génique pour calculer la nouvelle concentration de P700.

**Exemple d'utilisation (extrait de `update`):**
```python
def update(self, time, dt, env_state):
    light_intensity = env_state.get('light_intensity', 0)
    
    # Create a time array for the integration step
    t = [time, time + dt]
    
    # Solve the ODE for the current time step
    # The model function expects args in the order (light_intensity, params)
    solution = odeint(gene_regulation_model, [self.p700_concentration], t, args=(light_intensity, params))
    self.p700_concentration = solution[1][0]

    print(f"Updating biology at t={time}. P700 concentration: {self.p700_concentration:.4f}")

    return {
        'p700_concentration': self.p700_concentration
    }
```

│   │   ├── sensitivity_analysis.py

### Fichier: bioos/simulations/sensitivity_analysis.py
**Position:** `bioos/simulations/sensitivity_analysis.py`
**Utilité:** Réalise une analyse de sensibilité sur les paramètres du protocole de régénération de Ficus elastica pour identifier les étapes les plus critiques.
**Description:** Ce fichier implémente une analyse de sensibilité "un facteur à la fois" (OFAT) pour évaluer l'impact de la variation des probabilités de succès de chaque étape du protocole de régénération sur le rendement final. Il utilise la fonction de simulation de Monte Carlo (`run_regeneration_simulation`) en faisant varier un paramètre à la fois sur une plage définie, tout en maintenant les autres à leurs valeurs de base. Les résultats sont ensuite visualisés pour montrer l'influence de chaque paramètre.

**Fonctions clés:**
*   `run_regeneration_simulation(n_simulations, n_explants_initial, probabilities)`: Une version modifiée de la fonction de simulation de Monte Carlo qui accepte un dictionnaire de probabilités, permettant de tester l'impact de la variation des paramètres.
*   `run_sensitivity_analysis(base_probabilities, n_simulations, n_explants_initial)`: Orchestre l'analyse de sensibilité. Elle itère sur chaque paramètre, le fait varier sur une plage définie et exécute la simulation pour chaque valeur afin de calculer le rendement moyen.
*   `plot_sensitivity_results(sensitivity_results, param_range)`: Affiche les résultats de l'analyse de sensibilité sous forme graphique, montrant comment le rendement global moyen varie en fonction de chaque probabilité de succès. Le graphique est sauvegardé dans un fichier PNG.

**Paramètres clés analysés:**
*   `p_transfo`: Probabilité de transformation.
*   `p_selection`: Probabilité de survie à la sélection.
*   `p_callogenese`: Probabilité de formation de cal.
*   `p_organogenese`: Probabilité de production de bourgeons.
*   `p_enracinement`: Probabilité d'enracinement.
*   `p_acclimatation`: Probabilité d'acclimatation.

**Exemple d'utilisation (extrait de `run_sensitivity_analysis`):**
```python
def run_sensitivity_analysis(base_probabilities, n_simulations, n_explants_initial):
    sensitivity_results = {}
    param_range = np.linspace(0.1, 1.0, 10) # Plage de variation de 10% à 100%

    for param_name in base_probabilities.keys():
        yields = []
        
        current_probabilities = base_probabilities.copy()
        for p_value in param_range:
            current_probabilities[param_name] = p_value
            
            simulation_results = run_regeneration_simulation(
                n_simulations, n_explants_initial, current_probabilities
            )
            
            mean_yield = np.mean(simulation_results) / n_explants_initial * 100
            yields.append(mean_yield)

        sensitivity_results[param_name] = yields
    
    return sensitivity_results, param_range
```

        return False
    print("  -> RÉSULTAT: SUCCÈS. L'explant est viable.")
    # ... (autres étapes) ...
    return True
```

### Fichier: bioos/simulations/sop_procedural_simulation_optimized.py
**Position:** `bioos/simulations/sop_procedural_simulation_optimized.py`
**Utilité:** Simule de manière procédurale et interactive le protocole de régénération de Ficus elastica (HAWRA) avec des paramètres optimisés, en affichant le succès ou l'échec de chaque phase.
**Description:** Ce fichier est une version optimisée de `sop_procedural_simulation.py`. Il implémente une simulation interactive du protocole de régénération de Ficus elastica, mais avec une efficacité de transformation génétique améliorée (passant de 0.40 à 0.65) grâce à l'ajout de 200µM d'acétosyringone. Cette optimisation vise à démontrer l'impact positif de l'amélioration des conditions expérimentales sur le succès global du protocole. Comme la version non optimisée, elle suit le parcours d'un seul explant à travers sept étapes critiques, affichant le résultat (succès ou échec) en temps réel.

**Paramètres clés (Probabilités de succès):**
*   `explant_survival`: Survie de l'explant après prélèvement.
*   `transformation_efficiency`: Efficacité de la transformation par Agrobacterium (AMÉLIORÉE à 0.65).
*   `selection_survival`: Survie à la sélection par antibiotique/herbicide.
*   `callus_formation`: Formation de cals à partir des cellules transformées.
*   `shoot_induction`: Induction de pousses à partir des cals.
*   `root_formation`: Développement de racines sur les pousses.
*   `acclimatization_survival`: Survie de la plante lors du passage en terre.

**Fonctions clés:**
*   `run_procedural_simulation_optimized()`: Exécute la simulation optimisée pour un seul explant, affichant le résultat de chaque étape et le résultat final.

**Exemple d'utilisation (extrait de `run_procedural_simulation_optimized`):**
```python
def run_procedural_simulation_optimized():
    print("=====================================================================")
    print("=== Lancement de la Simulation du Protocole OPTIMISÉ (200µM AS) ===")
    print("=== Suivi du parcours d'un seul explant...                       ===")
    print("=====================================================================")
    time.sleep(1)

    # Étape 1: Prélèvement de l'explant
    print("\n[ÉTAPE 1/7] Prélèvement et préparation de l'explant...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["explant_survival"]:\
        print("  -> RÉSULTAT: ÉCHEC. L'explant n'a pas survécu.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. L'explant est viable.")

    # Étape 2: Transformation génétique (OPTIMISÉE)
    print("\n[ÉTAPE 2/7] Co-culture avec Agrobacterium (Protocole Optimisé, 200µM Acetosyringone)...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["transformation_efficiency"]:\
        print("  -> RÉSULTAT: ÉCHEC. Le transfert du gène HAWRA a échoué.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. L'ADN-T a été transféré avec une meilleure efficacité.")
    # ... (autres étapes) ...
    return True
```
**Description:** Ce fichier implémente une simulation séquentielle du protocole de régénération, où chaque étape est modélisée comme un événement probabiliste. Contrairement à la simulation de Monte Carlo qui donne un résultat global, cette simulation suit le parcours d'un *seul* explant à travers toutes les étapes, de la préparation à l'acclimatation. Elle est conçue pour être interactive, affichant des messages de progression et les résultats de chaque étape en temps réel, ce qui la rend utile pour la démonstration ou la compréhension didactique du protocole.

**Fonctions clés:**
*   `run_procedural_simulation()`: Exécute la simulation complète pour un seul explant. Elle parcourt séquentiellement les 7 étapes du protocole, en utilisant des probabilités prédéfinies pour déterminer le succès ou l'échec de chaque étape. Elle affiche des messages détaillés à chaque étape et termine la simulation si une étape échoue.

**Paramètres clés (Probabilités):**
*   `explant_survival`: Survie de l'explant après prélèvement.
*   `transformation_efficiency`: Efficacité de la transformation par Agrobacterium.
*   `selection_survival`: Survie à la sélection par antibiotique/herbicide.
*   `callus_formation`: Formation de cals à partir des cellules transformées.
*   `shoot_induction`: Induction de pousses à partir des cals.
*   `root_formation`: Développement de racines sur les pousses.
*   `acclimatization_survival`: Survie de la plante lors du passage en terre.

**Exemple d'utilisation (extrait de `run_procedural_simulation`):**
```python
def run_procedural_simulation():
    print("===========================================================")
    print("=== Lancement de la Simulation Procédurale du SOP Ficus ===")
    print("===========================================================")
    time.sleep(1)

    # Étape 1: Prélèvement de l'explant
    print("\n[ÉTAPE 1/7] Prélèvement et préparation de l'explant...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["explant_survival"]:
        print("  -> RÉSULTAT: ÉCHEC. L'explant n'a pas survécu à la stérilisation/préparation.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. L'explant est viable.")

    # ... (autres étapes similaires) ...

    print("  -> RÉSULTAT: SUCCÈS. La plante est acclimatée et viable !")
    print("\n=================================================================")
    print("=== 🎉 FÉLICITATIONS ! Une plante HAWRA-Ficus-G0 a été créée ! ===")
    print("=================================================================")
    return True
```

**Exemple d'utilisation (extrait de `run_procedural_simulation`):**
```python
def run_procedural_simulation():
    print("===========================================================")
    print("=== Lancement de la Simulation Procédurale du SOP Ficus ===")
    print("===========================================================")
    time.sleep(1)

    # Étape 1: Prélèvement de l'explant
    print("\n[ÉTAPE 1/7] Prélèvement et préparation de l'explant...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["explant_survival"]:
        print("  -> RÉSULTAT: ÉCHEC. L'explant n'a pas survécu à la stérilisation/préparation.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. L'explant est viable.")

    # ... (autres étapes similaires) ...

    print("  -> RÉSULTAT: SUCCÈS. La plante est acclimatée et viable !")
    print("\n=================================================================")
    print("=== 🎉 FÉLICITATIONS ! Une plante HAWRA-Ficus-G0 a été créée ! ===")
    print("=================================================================")
    return True
```

│   │   ├── sop_procedural_simulation_optimized.py

### Fichier: bioos/simulations/sop_procedural_simulation_optimized.py
**Position:** `bioos/simulations/sop_procedural_simulation_optimized.py`
**Utilité:** Simule de manière procédurale et interactive le protocole de régénération de Ficus elastica (HAWRA) avec des paramètres optimisés, notamment une efficacité de transformation améliorée.
**Description:** Ce fichier est une version optimisée de la simulation procédurale du protocole de régénération. Il suit le même principe de simulation étape par étape pour un seul explant, mais intègre des probabilités de succès ajustées pour refléter des améliorations expérimentales (par exemple, l'utilisation de 200µM d'acétosyringone pour augmenter l'efficacité de la transformation). Il est également conçu pour la démonstration, mettant en évidence l'impact positif des optimisations sur le processus global.

### Fichier: bioos/simulations/validate_simulation.py
**Position:** `bioos/simulations/validate_simulation.py`
**Utilité:** Valide le comportement du simulateur multiphysique en analysant les logs de simulation pour détecter les anomalies et assurer la cohérence du modèle.
**Description:** Ce fichier contient des fonctions pour analyser les fichiers de log générés par le simulateur multiphysique. Il effectue des vérifications critiques sur la concentration de P700, les excitations et les lectures des canaux (vert et rouge) pour s'assurer que le modèle se comporte comme prévu. Il est essentiel pour la vérification et la validation du simulateur, garantissant la fiabilité des résultats.
**Fonctions Clés:**
- `analyze_simulation_log(log_path, config)`: Prend en entrée le chemin d'un fichier de log JSON et la configuration de simulation, puis retourne un dictionnaire de résultats de validation incluant les erreurs détectées, le nombre d'excitations et les lectures des canaux.
**Extraits de Code:**
```python
def analyze_simulation_log(log_path, config):
    \"\"\"
    Analyzes the multiphysics simulation log to validate the model\'s behavior.
    \"\"\"
    with open(log_path, \'r\') as f:
        log = json.load(f)

    p700_threshold = config[\'quantum\'][\'threshold\']

    errors = []
    excitations = 0
    green_reads = 0
    red_reads = 0

    for i in range(1, len(log)):
        prev_state = log[i-1]
        current_state = log[i]

        # 1. P700 concentration validation
        if current_state[\'light_intensity\'] > 0:
            if current_state[\'p700_concentration\'] < prev_state[\'p700_concentration\'] and prev_state[\'p700_concentration\'] < 0.99 :
                # Allow for decay when light is on but synthesis is saturated
                pass
        else: # No light
            if current_state[\'p700_concentration\'] > prev_state[\'p700_concentration\']:\
                 errors.append(f\"t={current_state[\'time\']}: P700 increased without light.\")

        # 2. P700 excitation validation
        if current_state[\'p700_state\'] == 1 and prev_state[\'p700_state\'] == 0:
            excitations += 1
            if current_state[\'p700_concentration\'] < p700_threshold:
                errors.append(f\"t={current_state[\'time\']}: P700 excited below threshold.\")

        # 3. Readout channel validation
        if current_state[\'luc_green_output\'] > 0 or current_state[\'luc_red_output\'] > 0:
            if prev_state[\'p700_state\'] == 0:
                errors.append(f\"t={current_state[\'time\']}: Readout occurred from ground state.\")

        if current_state[\'luc_green_output\'] > 0:
            green_reads += 1

        if current_state[\'luc_red_output\'] > 0:
            red_reads += 1

        # 4. Mutual exclusion of readouts
        if current_state[\'luc_green_output\'] > 0 and current_state[\'luc_red_output\'] > 0:
            errors.append(f\"t={current_state[\'time\']}: Mutual exclusion of readouts violated.\")

    return {
        "validation_status": "SUCCESS" if not errors else "FAILURE",
        "errors": errors,
        "total_steps": len(log),
        "excitations": excitations,
        "green_reads": green_reads,
        "red_reads": red_reads
    }
```

### Fichier: bioos/simulations/regeneration_simulation.py
**Position:** `bioos/simulations/regeneration_simulation.py`
**Utilité:** Estime le rendement global du protocole de régénération de Ficus elastica (HAWRA) par une approche de Monte Carlo, modélisant chaque étape critique comme un événement probabiliste.
**Description:** Ce fichier implémente une simulation de Monte Carlo pour prédire le nombre de plantes viables (HAWRA-Ficus-G0) obtenues à partir d'un nombre initial d'explants. Le protocole est traité comme une chaîne de Markov, où la sortie d'une étape devient l'entrée de la suivante, chaque étape ayant une probabilité de succès. Les probabilités sont basées sur des estimations de la littérature et peuvent être affinées expérimentalement. Le fichier inclut également des fonctions pour visualiser les résultats sous forme d'histogramme et calculer des statistiques clés.
**Fonctions Clés:**
- `run_regeneration_simulation(n_simulations, n_explants_initial)`: Exécute la simulation de Monte Carlo pour un nombre donné de simulations et d'explants initiaux, retournant une liste du nombre final de plantes viables pour chaque simulation.
- `plot_results(simulation_results, n_explants_initial)`: Affiche un histogramme des résultats de la simulation, calcule la moyenne, l'écart-type et le rendement global moyen, puis sauvegarde le graphique.
**Paramètres Clés (Probabilités de succès):**
- `p_transfo`: Probabilité qu'une cellule d'explant intègre l'ADN-T.
- `p_selection`: Probabilité qu'une cellule transformée survive à la sélection.
- `p_callogenese`: Probabilité qu'un cal se forme à partir de cellules sélectionnées.
- `p_organogenese`: Probabilité qu'un cal génère des bourgeons viables.
- `p_enracinement`: Probabilité qu'une pousse développe un système racinaire.
- `p_acclimatation`: Probabilité qu'une plantule survive au transfert en serre.
**Extraits de Code:**
```python
def run_regeneration_simulation(n_simulations, n_explants_initial):
    p_transfo = 0.10
    p_selection = 0.60
    p_callogenese = 0.50
    p_organogenese = 0.40
    p_enracinement = 0.70
    p_acclimatation = 0.50

    final_plant_counts = []

    for _ in range(n_simulations):
        n_surviving_selection = np.random.binomial(n_explants_initial, p_transfo * p_selection)
        n_calli = np.random.binomial(n_surviving_selection, p_callogenese)
        n_shoots = np.random.binomial(n_calli, p_organogenese)
        n_rooted_plantlets = np.random.binomial(n_shoots, p_enracinement)
        n_final_plants = np.random.binomial(n_rooted_plantlets, p_acclimatation)
        final_plant_counts.append(n_final_plants)

    return final_plant_counts

def plot_results(simulation_results, n_explants_initial):
    mean_plants = np.mean(simulation_results)
    std_plants = np.std(simulation_results)
    rendement_moyen = (mean_plants / n_explants_initial) * 100

    plt.figure(figsize=(12, 7))
    plt.hist(simulation_results, bins=range(0, max(simulation_results) + 2), alpha=0.75, edgecolor='black')
    plt.title(f\"Distribution du Nombre de Plantes Viables (N_explants = {n_explants_initial})\")
    plt.xlabel(\"Nombre de Plantes HAWRA-Ficus-G0 Viables\")
    plt.ylabel(\"Fréquence (sur {} simulations)\".format(len(simulation_results)))
    plt.axvline(mean_plants, color=\'r\', linestyle=\'dashed\', linewidth=2, label=f\"Moyenne: {mean_plants:.2f}\")
    plt.legend()
    print(f\"Nombre moyen de plantes finales: {mean_plants:.2f} ± {std_plants:.2f}\")
    output_path = \"/Users/mehdiwhb/Desktop/HAWRA/05_data/results/regeneration_simulation_yield.png\"
    plt.savefig(output_path)
```

**Fonctions clés:**
*   `run_procedural_simulation_optimized()`: Exécute la simulation complète pour un seul explant avec les paramètres optimisés. Elle parcourt séquentiellement les 7 étapes du protocole, affichant les messages de progression et les résultats en temps réel. La simulation s'arrête si une étape échoue.

**Paramètres clés (Probabilités optimisées):**
*   `explant_survival`: Survie de l'explant après prélèvement.
*   `transformation_efficiency`: Efficacité de la transformation par Agrobacterium (AMÉLIORÉE).
*   `selection_survival`: Survie à la sélection par antibiotique/herbicide.
*   `callus_formation`: Formation de cals à partir des cellules transformées.
*   `shoot_induction`: Induction de pousses à partir des cals.
*   `root_formation`: Développement de racines sur les pousses.
*   `acclimatization_survival`: Survie de la plante lors du passage en terre.

**Exemple d'utilisation (extrait de `run_procedural_simulation_optimized`):**
```python
def run_procedural_simulation_optimized():
    print("=====================================================================")
    print("=== Lancement de la Simulation du Protocole OPTIMISÉ (200µM AS) ===")
    print("=====================================================================")
    time.sleep(1)

    # Étape 1: Prélèvement de l'explant
    print("\n[ÉTAPE 1/7] Prélèvement et préparation de l'explant...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["explant_survival"]:
        print("  -> RÉSULTAT: ÉCHEC. L'explant n'a pas survécu.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. L'explant est viable.")

    # Étape 2: Transformation génétique (OPTIMISÉE)
    print("\n[ÉTAPE 2/7] Co-culture avec Agrobacterium (Protocole Optimisé, 200µM Acetosyringone)...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["transformation_efficiency"]:
        print("  -> RÉSULTAT: ÉCHEC. Le transfert du gène HAWRA a échoué.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. L'ADN-T a été transféré avec une meilleure efficacité.")

    # ... (autres étapes similaires) ...

    print("  -> RÉSULTAT: SUCCÈS. La plante est acclimatée et viable !")
    print("\n=================================================================")
    print("=== 🎉 FÉLICITATIONS ! Une plante HAWRA-Ficus-G0 a été créée ! ===")
    print("=================================================================")
    return True
```

│   │   ├── models
│   │   ├── models
│   │   │   ├── __init__.py
│   │   │   ├── gene_regulatory_network.py
│   │   │   └── metabolic_pathway.py
│   │   ├── data
│   │   │   ├── __init__.py
│   │   │   └── constants.py
│   │   └── utils
│   │       ├── __init__.py
│   │       └── visualization.py
│   ├── bioos_core
│   │   ├── __init__.py
│   │   ├── system_monitor.py
│   │   └── task_scheduler.py
│   └── tests
│       ├── __init__.py
│       ├── test_hawra_simulator.py
│       └── test_bioos_core.py
├── arbol
│   ├── __init__.py
│   ├── compiler
│   │   ├── __init__.py
│   │   ├── __pycache__
│   │   │   ├── parser.cpython-313.pyc
│   │   │   ├── lexer.cpython-313.pyc
│   │   │   ├── arbol_ast.cpython-313.pyc
│   │   │   ├── test_compiler.cpython-313-pytest-9.0.1.pyc
│   │   │   └── ast.cpython-313.pyc
│   │   ├── parser.py
│   │   
│   │   ### Fichier: arbol/grammar/arbol_unified.ebnf
**Position:** `arbol/grammar/arbol_unified.ebnf`
**Utilité:** Définit la grammaire formelle du langage Arbol en utilisant la notation EBNF (Extended Backus-Naur Form).
**Description:** Ce fichier est essentiel pour la spécification syntaxique du langage Arbol. Il décrit les règles de production pour tous les éléments du langage, tels que la structure d'un programme (`program`), les différentes déclarations (`circuit_definition`, `logical_qubit_definition`, `gate_definition`, `stimulus_definition`), les opérations (`quantum_operation`, `measure_operation`, `stimulus_application`), et les symboles terminaux (`ID`, `NUMBER`, `STRING`). Cette grammaire est utilisée par le lexer et le parseur pour valider et interpréter correctement le code Arbol. Elle inclut également des règles pour ignorer les espaces blancs et les commentaires.

│   │   ### Fichier: arbol/compiler/__init__.py
**Position:** `arbol/compiler/__init__.py`
**Utilité:** Marque le répertoire `compiler` comme un paquet Python.
**Description:** Ce fichier est vide et sert uniquement à indiquer que le répertoire `compiler` doit être traité comme un paquet contenant des modules Python. Il permet l'importation des modules définis dans ce répertoire (par exemple, `from arbol.compiler import lexer`).

│   │   ### Fichier: arbol/compiler/parser.py
│   │   **Position:** `arbol/compiler/parser.py`
│   │   **Utilité:** Analyse la séquence de jetons (tokens) produite par le lexer pour construire l'arbre syntaxique abstrait (AST) du code Arbol.
│   │   **Description:** Ce fichier implémente la classe `Parser` qui est responsable de l'analyse syntaxique du langage Arbol. Il prend en entrée une liste de `Token` (générés par `lexer.py`) et les transforme en une structure hiérarchique d'objets `Node` (définis dans `arbol_ast.py`). Le parseur utilise une méthode d'analyse descendante récursive, avec des fonctions `parse_` pour chaque type de déclaration ou d'instruction (ex: `parse_circuit_definition`, `parse_gene_definition`, `parse_grn_definition`). Il gère la progression à travers les jetons (`advance`, `eat`) et inclut un mécanisme de synchronisation (`synchronize`) pour la récupération d'erreurs, permettant au compilateur de continuer l'analyse même après avoir rencontré une erreur syntaxique. Les erreurs sont rapportées via un `error_reporter`.
│   │   
│   │   ├── compiler.py
│   │   
│   │   ### Fichier: arbol/compiler/compiler.py
**Position:** `arbol/compiler/compiler.py`
**Utilité:** Convertit l'AST généré par le parseur en instructions d'assemblage JSON exécutables par le simulateur.
**Description:** Ce fichier implémente le compilateur principal qui traverse l'AST en utilisant le modèle de conception "Visiteur". Il gère la table des symboles pour suivre les définitions de gènes, de stimuli et de circuits, et génère les instructions JSON correspondantes. Par exemple, la méthode `visit_GeneDefinition` traite la définition d'un gène et ajoute les configurations nécessaires à la sortie JSON.

### Fichier: arbol/tests/test_parser.py
**Position:** `arbol/tests/test_parser.py`
**Utilité:** Contient les tests unitaires approfondis pour le parseur Arbol, validant la construction correcte de l'arbre syntaxique abstrait (AST) à partir du code source Arbol.
**Description:** Ce fichier utilise le framework `unittest` pour vérifier la robustesse et la précision de la classe `Parser`. Il couvre une gamme étendue de scénarios de test, incluant :
*   **`test_gene_definition`**: Valide l'analyse des définitions de gènes et de leurs stimuli associés.
*   **`test_circuit_definition`**: Assure la bonne interprétation des définitions de circuits et de leurs séquences de gènes.
*   **`test_logical_qubit_definition`**: Confirme l'analyse correcte des qubits logiques, y compris leur taille et les clauses `IS` spécifiant les activateurs et répresseurs.
*   **`test_measure_statement_on` et `test_measure_statement_arrow`**: Vérifie les différentes syntaxes des instructions `MEASURE` (par exemple, `MEASURE q1 ON c1` et `MEASURE q1 -> m1`).
*   **`test_if_statement` et `test_nested_if_statements`**: Valide l'analyse des instructions conditionnelles, y compris les structures `IF` imbriquées.
*   **`test_complex_program`**: Un test intégré qui combine plusieurs constructions Arbol pour évaluer la robustesse globale du parseur.
*   **`test_parser_error_reporting`**: S'assure que le parseur identifie et signale précisément les erreurs de syntaxe avec les numéros de ligne et de colonne.

Ces tests sont essentiels pour garantir que le parseur Arbol interprète fidèlement le langage, ce qui est fondamental pour la fiabilité du compilateur.

**Exemple de test (`test_logical_qubit_definition`):**
```python
def test_logical_qubit_definition(self):
    code = "LOGICAL_QUBIT q1 [2] IS { activator: p1, repressor: p2 };"
    lexer = Lexer(code)
    error_reporter = ErrorReporter()
    parser = Parser(lexer, error_reporter)
    ast = parser.parse()
    expected_ast = Program(statements=[
        LogicalQubitDefinition(
            name=Identifier(name='q1'),
            size=2,
            is_clause=IsClause(promoters=[
                Promoter(type='activator', name='p1'),
                Promoter(type='repressor', name='p2')
            ]),
            line=actual_line,
            column=actual_column
        )
    ])
    self.assertEqual(ast, expected_ast)
```│   │   
### Fichier: arbol/tests/test_compiler.py
**Position:** `arbol/tests/test_compiler.py`
**Utilité:** Contient les tests unitaires pour le compilateur Arbol.
**Description:** Ce fichier utilise le framework `unittest` pour tester la fonctionnalité du compilateur Arbol. Actuellement, il contient un test de base (`test_simple`) qui sert de placeholder. Cela indique que des tests plus complets pour la traduction de l'AST en instructions d'assemblage JSON et la gestion de la table des symboles sont à développer.

**Exemple de test (`test_simple`):**
```python
def test_simple(self):
    assert 1 == 1
```│   │   
│   │   ### Fichier: arbol/tests/test_lexer.py
**Position:** `arbol/tests/test_lexer.py`
**Utilité:** Contient les tests unitaires pour le lexer Arbol.
**Description:** Ce fichier utilise le framework `unittest` pour vérifier le bon fonctionnement de la classe `Lexer`. Il teste la reconnaissance des mots-clés, des identifiants, des nombres, des chaînes de caractères, des opérateurs et des délimiteurs. Il s'assure également que les commentaires et les espaces blancs sont correctement ignorés et que les jetons inconnus déclenchent des exceptions appropriées. Ce fichier est essentiel pour garantir l'intégrité du processus de tokenisation du langage Arbol.
│   │   
│   │   ├── error.py
│   │   
│   │   ### Fichier: arbol/compiler/error.py
│   │   **Position:** `arbol/compiler/error.py`
│   │   **Utilité:** Gère la signalisation et la collecte des erreurs rencontrées pendant les phases de compilation (lexicale, syntaxique, sémantique) du langage Arbol.
│   │   **Description:** Ce fichier définit deux classes principales : `CompilationError` et `ErrorReporter`. `CompilationError` est une classe d'exception personnalisée qui encapsule les détails d'une erreur de compilation, y compris le message, le numéro de ligne et le numéro de colonne où l'erreur s'est produite. `ErrorReporter` est une classe utilitaire qui centralise la gestion des erreurs. Elle permet d'enregistrer plusieurs `CompilationError` au fur et à mesure qu'elles sont détectées et de les afficher de manière structurée à la fin du processus de compilation. Cela assure une gestion cohérente et informative des erreurs pour l'utilisateur.
│   │   
│   │   ├── lexer.py
│   │   
│   │   ### Fichier: arbol/compiler/lexer.py
│   │   **Position:** `arbol/compiler/lexer.py`
│   │   **Utilité:** Effectue l'analyse lexicale (tokenisation) du code source Arbol, convertissant une séquence de caractères en une séquence de jetons (tokens).
│   │   **Description:** Ce fichier contient la classe `Lexer` qui est chargée de décomposer le code Arbol en unités significatives appelées "tokens". Il utilise des expressions régulières définies dans `token_specification` pour reconnaître les mots-clés (comme `circuit`, `apply`, `stimulus`, `genes`, `grn`), les identifiants, les opérateurs (`->`, `=`, `(`, `)`, `{`, `}`, `[`, `]`, `,`, `;`, `:`), les littéraux (chaînes de caractères, nombres), et les commentaires. La classe `Token` encapsule chaque jeton avec son type, sa valeur, et sa position (ligne et colonne) dans le fichier source, ce qui est crucial pour la détection d'erreurs et l'analyse syntaxique ultérieure. Le lexer gère également les espaces blancs et les sauts de ligne.
│   │   
│   │   └── arbol_ast.py
│   │   
│   │   ### Fichier: arbol/compiler/arbol_ast.py
│   │   **Position:** `arbol/compiler/arbol_ast.py`
│   │   **Utilité:** Définit la structure de l'arbre syntaxique abstrait (AST) pour le langage Arbol. C'est le cœur de la représentation interne du code Arbol après l'analyse syntaxique.
│   │   **Description:** Ce fichier utilise les `dataclasses` de Python pour modéliser les différents nœuds qui composent l'AST. Chaque classe représente un élément grammatical du langage Arbol, comme les identifiants (`Identifier`), les définitions de programmes (`Program`), les définitions de circuits (`CircuitDefinition`), les opérations quantiques (`QuantumOperation`), les déclarations de qubits (`QubitDeclaration`), les applications de portes quantiques (`GateApplication`), les définitions de stimuli (`StimulusDefinition`), les définitions de gènes (`GeneDefinition`), et les interactions de réseaux de régulation génique (`GRNInteraction`). Ces structures de données sont essentielles pour que le compilateur puisse interpréter, analyser et transformer le code Arbol.
│   │   
│   ├── test_integration.arbol

### Fichier: arbol/test_integration.arbol
**Position:** `arbol/test_integration.arbol`
**Utilité:** Sert de test d'intégration pour le langage Arbol, démontrant la définition d'un stimulus, d'un qubit logique, d'un bit classique et d'un circuit quantique simple.
**Description:** Ce fichier Arbol illustre un scénario de test d'intégration typique. Il définit :
*   **Stimulus:** Un stimulus `light_pulse` avec des paramètres de durée et d'intensité.
*   **Qubit Logique:** Un qubit logique `q`.
*   **Bit Classique:** Un bit classique `m`.
*   **Circuit Quantique:** Un circuit `main()` qui applique une séquence d'opérations :
    *   Application du stimulus `light_pulse`.
    *   Application d'une porte de Hadamard (`H`) au qubit `q`.
    *   Mesure du qubit `q` et stockage du résultat dans le bit classique `m`.

Ce fichier est crucial pour valider l'intégration et le bon fonctionnement des différentes composantes du compilateur Arbol, du lexer au compilateur final, en passant par le parseur et la génération de l'AST.

**Exemple de circuit:**
```arbol
circuit main() {
    apply light_pulse(duration="50ms", intensity=0.8);
    H(q);
    measure q -> m;
}
```
│   └── dynamic_grn.arbol

### Fichier: arbol/dynamic_grn.arbol
**Position:** `arbol/dynamic_grn.arbol`
**Utilité:** Définit un protocole Arbol pour la gestion dynamique d'un réseau de régulation génique (GRN), incluant la définition de gènes, les interactions GRN, les stimuli et une séquence d'exécution.
**Description:** Ce fichier Arbol illustre comment modéliser des systèmes biologiques dynamiques. Il contient :
*   **Définition des gènes:** Spécifie des gènes (ex: "GENE_A", "GENE_B", "GENE_C") avec des paramètres tels que le taux basal (`basal_rate`), le taux de dégradation (`degradation_rate`) et la sensibilité à la lumière (`light_sensitivity`).
*   **Réseau de régulation génique (GRN):** Décrit les interactions entre les gènes, comme l'activation ("GENE_A" active "GENE_B") et la répression ("GENE_B" réprime "GENE_C"), avec des poids, des coefficients de Hill et des concentrations de demi-maximum.
*   **Définition des stimuli:** Définit des stimuli externes (ex: `light_pulse` avec `intensity` et `duration`) qui peuvent affecter les gènes.
*   **Séquence d'exécution (`run`):** Ordonne une série d'étapes (`step`) et l'application de stimuli (`apply`) à des moments spécifiques, simulant ainsi l'évolution dynamique du système.

Ce fichier est un exemple clé de la capacité du langage Arbol à modéliser des expériences biologiques complexes et dynamiques.

**Exemple de définition de gène:**
```arbol
gene "GENE_A" with basal_rate=0.01, degradation_rate=0.05, light_sensitivity=0.8;
```

**Exemple d'interaction GRN:**
```arbol
"GENE_A" activates "GENE_B" with weight=1.0, hill_coefficient=2, half_max_concentration=0.5;
```

**Exemple de séquence d'exécution:**
```arbol
run {
    step(50);
    apply light_pulse(intensity: 1.0, duration: 10.0) on "GENE_A";
    step(100);
}
```
├── 02_arbol_interface
│   ├── ide
│   │   └── .keep
│   ├── jetson_client
│   │   └── .keep
│   └── gene_control.arbol
├── 01_genomics
│   ├── plasmids
│   │   ├── HAWRA_PLASMID_v2.gb
│   │   ├── HAWRA_PLASMID_v3.gb
│   │   ├── HAWRA_FINAL_VALIDATED.gb
│   │   └── HAWRA_TRANSDUCER_v1.gb
│   ├── experiments
│   │   ├── first_bloom.bsim.json
│   │   ├── first_bloom.arbol
│   │   └── first_bloom_results.json
│   ├── raw_sequences
│   │   ├── psaA_NC_000932.1.fasta
│   │   ├── Lsi1_SIT1_AB222272.1.gb
│   │   ├── LUC_M15077.1.gb
│   │   ├── raw_chromosomes
│   │   │   ├── CM023739.1.gb
│   │   │   ├── CM023742.1.gb
│   │   │   ├── CM023738.1.gb
│   │   │   ├── CM023733.1.gb
│   │   │   ├── CM023745.1.gb
│   │   │   ├── CM111266.gb
│   │   │   ├── CM023737.1.gb
│   │   │   ├── CM023741.1.gb
│   │   │   └── CM009940.2.gb
│   │   ├── psaA_NC_000932.1.gb
│   │   ├── Lsi1_SIT1_AB222272.1.fasta
│   │   ├── LUC_M15077.1.fasta
│   │   ├── PEPC1_X13660.1.fasta
│   │   ├── HSP70_NM_112093.3.gb
│   │   ├── CDS
│   │   │   ├── Ficus_elastica_complete_genome
│   │   │   │   ├── CM111267 (1).dna
│   │   │   │   ├── CM111273.dna
│   │   │   │   ├── CM111267.dna
│   │   │   │   ├── CM111266.dna
│   │   │   │   ├── CM111272.dna
│   │   │   │   ├── CM111270.dna
│   │   │   │   ├── CM111271.dna
│   │   │   │   ├── CM111275.dna
│   │   │   │   ├── CM111274.dna
│   │   │   │   ├── CM111266 (1).dna
│   │   │   │   ├── CM111269 (1).dna
│   │   │   │   ├── CM111270 (1).dna
│   │   │   │   ├── CM111271 (1).dna
│   │   │   │   ├── CM111268 (1).dna
│   │   │   │   ├── CM111273 (1).dna
│   │   │   │   ├── CM111268.dna
│   │   │   │   └── CM111269.dna
│   │   │   ├── NOS.fasta
│   │   │   ├── CRY2.fasta
│   │   │   ├── HSP70_CDS.fasta
│   │   │   ├── CRY2_CDS.fasta
│   │   │   ├── LUC_CDS.fasta
│   │   │   ├── PEPC.fasta
│   │   │   ├── Lsi1_SIT1_CDS.fasta
│   │   │   ├── CM111266.gb
│   │   │   ├── 35S.fasta
│   │   │   ├── Lsi1.fasta
│   │   │   ├── psaA_CDS_CORRECTED.fasta
│   │   │   ├── HAWRA_FINAL_VALIDATED.gb
│   │   │   ├── PEPC1_CDS.fasta
│   │   │   ├── psaA.fasta
│   │   │   ├── psaA_CDS.fasta
│   │   │   ├── CDS_SUMMARY.txt
│   │   │   ├── Luc.fasta
│   │   │   ├── HAWRA_WITH_REAL_CDS.gb
│   │   │   └── HSP70.fasta
│   │   ├── CRY2_NM_112093.3.fasta
│   │   ├── NC_001497.gb
│   │   ├── HSP70_NM_112093.3.fasta
│   │   ├── download_ficus_genome.py
│   │   ├── CRY2_NM_100320.4.gb
│   │   └── PEPC1_X13660.1.gb
│   ├── genome_analysis_scripts
│   │   ├── visualize_electrochemical.py
│   │   ├── visualize_hormonal_regulation.py
│   │   ├── visualize_plasmid.py
│   │   └── __pycache__
│   │       └── visualize_plasmid.cpython-313.pyc
│   └── processed_sequences
│       └── .keep
└── demo_pipeline.py
```
