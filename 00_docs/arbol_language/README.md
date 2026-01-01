# HAWRA - ARBOL Language

## Introduction

ARBOL (A-life Regulation and Bio-computation Orchestration Language) est un langage de programmation conçu pour interagir avec des systèmes bio-hybrides, en particulier le projet HAWRA.

Contrairement aux langages de programmation traditionnels qui s'exécutent sur du silicium, ARBOL est conçu pour être compilé en signaux biologiques et en constructions génétiques qui opèrent au sein d'un organisme vivant.

## Philosophie

La philosophie d'ARBOL est basée sur les principes suivants :

- **Orchestration multi-échelle** : Le langage doit permettre de décrire des opérations à différents niveaux de l'organisme : moléculaire, cellulaire et tissulaire.
- **Computation incarnée** : Les programmes ARBOL ne sont pas des simulations, mais des instructions qui modifient directement l'état et le comportement de l'organisme.
- **Inspiration biologique** : La syntaxe et les concepts du langage s'inspirent des processus biologiques réels, tels que la régulation génique, la signalisation cellulaire et la morphogenèse.

Ce document sert de point de départ pour la formalisation du langage ARBOL.

## Sortie de compilation (BSIM)

Le compilateur Arbol génère un script de bio‑simulation au format BSIM (JSON): `*.bsim.json`.

- Guide détaillé: voir `BSIM_format.md` dans ce dossier.
- Contenu: `metadata`, liste `instructions` (`INITIALIZE`, `QUANTUM_OP`, `MEASURE`, `RUN_UNTIL`, `run`).
- Utilisation: `python3 -m 04_arbol.compiler.compiler 04_arbol/<script>.arbol`.

## Syntaxe de Base

ARBOL utilise une syntaxe simple et déclarative. Voici les éléments principaux :

*   **Déclaration de Qubits** : Vous pouvez allouer des qubits en spécifiant leur type et leur nombre.
    ```arbol
    // Alloue 10 qubits de type P700
    qubits = alloc(P700, 10)
    ```

*   **Application de Portes** : Appliquez des portes quantiques aux qubits alloués.
    ```arbol
    // Applique une porte Hadamard au premier qubit
    H(qubits[0])

    // Applique une porte CNOT entre le premier et le deuxième qubit
    CNOT(qubits[0], qubits[1])
    ```

*   **Mesure** : Mesurez l'état des qubits.
    ```arbol
    // Mesure tous les qubits
    result = measure(qubits)
    ```

## Exemple Complet : Algorithme de Grover

Voici un exemple de l'algorithme de Grover pour rechercher un état spécifique dans une base de données de 4 éléments.

```arbol
// 1. Initialisation
qubits = alloc(P700, 2)
H(qubits[0])
H(qubits[1])

// 2. Oracle (marque l'état |11>)
CZ(qubits[0], qubits[1])

// 3. Amplification
H(qubits[0])
H(qubits[1])
X(qubits[0])
X(qubits[1])
CZ(qubits[0], qubits[1])
X(qubits[0])
X(qubits[1])
H(qubits[0])
H(qubits[1])

// 4. Mesure
result = measure(qubits)
```

Ce code est ensuite compilé par le compilateur ARBOL en une séquence d'instructions que le BioOS peut exécuter sur la PQPE.
