# Rapport d'Analyse de Sensibilité du Protocole de Régénération HAWRA

**Date:** 2024-07-26
**Auteur:** Gemini, Assistant de Codage IA
**Statut:** Final

## 1. Objectif

Ce rapport présente les résultats d'une analyse de sensibilité numérique menée sur le protocole de régénération de *Ficus elastica* pour le projet HAWRA. L'objectif était d'identifier et de hiérarchiser les étapes les plus critiques du protocole, afin de guider de manière stratégique les futurs efforts d'optimisation expérimentale.

## 2. Méthodologie

Nous avons employé une analyse de sensibilité "un facteur à la fois" (OFAT) basée sur une simulation de Monte Carlo du protocole. Le modèle simule le parcours d'une cohorte d'explants à travers six étapes probabilistes clés :

1.  **Transformation** (`p_transfo`)
2.  **Sélection** (`p_selection`)
3.  **Callogenèse** (`p_callogenese`)
4.  **Organogenèse** (`p_organogenese`)
5.  **Enracinement** (`p_enracinement`)
6.  **Acclimatation** (`p_acclimatation`)

Pour chaque paramètre, sa probabilité de succès a été variée de 10% à 100% (par pas de 10%), tandis que tous les autres paramètres étaient maintenus à leur valeur de base. Pour chaque point de données, 2000 simulations complètes ont été exécutées pour calculer un rendement global moyen robuste.

## 3. Résultats

L'analyse a généré un ensemble de données quantifiant la relation entre chaque paramètre et le rendement global. Les résultats sont visualisés dans le graphique ci-dessous.

![Analyse de Sensibilité]( /Users/mehdiwhb/Desktop/HAWRA/05_data/results/sensitivity_analysis.png)

*Figure 1: Graphique de l'analyse de sensibilité montrant l'impact de la variation de chaque probabilité d'étape sur le rendement global moyen du protocole.*

## 4. Discussion et Hiérarchisation des Facteurs Critiques

Le graphique révèle une hiérarchie claire et prononcée des facteurs limitants. La sensibilité du rendement global à chaque paramètre, indiquée par la pente de la courbe correspondante, est la suivante (par ordre décroissant d'influence) :

1.  **`p_transfo` (Efficacité de la transformation)** : Ce paramètre est, de manière écrasante, le plus critique. En raison de sa position en amont de la chaîne de processus et de sa faible valeur de base (10%), toute variation de `p_transfo` a un effet multiplicateur majeur sur le résultat final. Le rendement est presque directement proportionnel à l'efficacité de la transformation.

2.  **`p_organogenese` (Organogenèse)** : Le deuxième facteur le plus influent. L'amélioration de la capacité des cals à générer des bourgeons a un impact significatif sur le rendement.

3.  **`p_callogenese` (Callogenèse)** : Le troisième facteur critique. La formation de cals à partir des cellules sélectionnées est un prérequis essentiel pour les étapes ultérieures.

4.  **`p_acclimatation` (Acclimatation)** : La survie des plantules lors du transfert en serre est le quatrième facteur le plus sensible.

5.  **`p_selection` (Sélection)** : Bien qu'importante, l'efficacité de la sélection a un impact moindre que les étapes précédentes dans la plage étudiée.

6.  **`p_enracinement` (Enracinement)** : Ce paramètre est le moins sensible. Bien qu'un enracinement réussi soit nécessaire, des améliorations dans ce domaine produiront des gains de rendement marginaux par rapport aux autres facteurs.

## 5. Recommandations Stratégiques

Sur la base de cette validation numérique, nous formulons les recommandations suivantes pour la prochaine phase d'optimisation expérimentale :

*   **Priorité Absolue :** Les efforts de R&D doivent être massivement concentrés sur l'amélioration de l'efficacité de la transformation (`p_transfo`). Les stratégies pourraient inclure :
    *   L'optimisation de la concentration d'acetosyringone (comme déjà testé).
    *   L'exploration de différentes souches d'*Agrobacterium*.
    *   Le test de méthodes de transformation alternatives (biolistique, etc.).

*   **Priorité Secondaire :** Une fois des gains significatifs obtenus sur la transformation, les efforts devraient se tourner vers l'optimisation de l'organogenèse et de la callogenèse, par exemple en ajustant les régulateurs de croissance dans les milieux de culture.

*   **Faible Priorité :** L'optimisation de l'enracinement ne devrait être entreprise qu'une fois que les goulots d'étranglement majeurs auront été résolus.

Cette approche hiérarchisée, guidée par les données, garantit une allocation efficace des ressources et maximise les chances d'améliorer de manière significative le rendement global du protocole HAWRA.
