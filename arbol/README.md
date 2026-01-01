# Arbol: Un compilateur pour la biologie synthétique quantique

Arbol est un compilateur pour un langage de haut niveau qui décrit des circuits génétiques et des opérations quantiques. Il produit un script de bio-simulation au format JSON qui peut être exécuté par une plateforme de simulation.

## Langage Arbol

Le langage Arbol permet de définir :

- **Stimuli**: Des signaux externes qui peuvent influencer le système (ex: lumière, chaleur).
- **Qubits logiques**: Des qubits définis par des constructions génétiques.
- **Portes quantiques**: Des opérations quantiques personnalisées.
- **Bits classiques**: Des bits classiques pour stocker les résultats de mesure.
- **Circuits**: Des séquences d'opérations quantiques.
- **Opérations de mesure**: Des mesures de qubits qui stockent le résultat dans des bits classiques.

### Exemple de code Arbol

```arbol
(* Définition d'un stimulus lumineux *)
stimulus light_pulse is light(intensity: float, duration: string);

(* Définition d'un qubit logique *)
logical_qubit q1 is {
    activator: geneA,
    repressor: geneB
};

(* Définition d'un bit classique *)
classical_bit c1;

(* Définition d'une porte quantique personnalisée *)
gate h q1 {
    H(q1)
}

(* Définition d'un circuit quantique *)
circuit my_circuit {
    q1, q2
    H(q1)
    CNOT(q1, q2)
}

(* Application d'un stimulus *)
apply light_pulse to q1 with {
    intensity = 1.0,
    duration = "10s"
}

(* Mesure d'un qubit *)
measure q1 to c1;
```

## Compilateur

Le compilateur Arbol est écrit en Python et se trouve dans le répertoire `compiler/`.

### Utilisation

Pour compiler un fichier Arbol (`.arbol`), utilisez la commande suivante :

```bash
python3 compiler/compiler.py <input_file.arbol>
```

Le compilateur générera un fichier de bio-simulation JSON (`.bsim.json`) dans le même répertoire.

### Format BSIM

- Sortie par défaut: `*.bsim.json` (contenu JSON).
- Structure: `metadata` + `instructions` (`INITIALIZE`, `QUANTUM_OP`, `MEASURE`, `RUN_UNTIL`, `run`).
- Renommage possible en `*.bsim` (le contenu reste du JSON).
- Documentation détaillée: `00_docs/arbol_language/BSIM_format.md`.

### Tests

Les tests du compilateur se trouvent dans `compiler/test_compiler.py`. Pour les exécuter, utilisez la commande suivante :

```bash
python3 -m unittest compiler/test_compiler.py
```
